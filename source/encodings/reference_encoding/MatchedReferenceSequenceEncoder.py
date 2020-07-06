import abc
import os

from source.IO.dataset_export.PickleExporter import PickleExporter
from source.caching.CacheHandler import CacheHandler
from source.encodings.DatasetEncoder import DatasetEncoder
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.SequenceMatchingSummaryType import SequenceMatchingSummaryType
from source.util.ParameterValidator import ParameterValidator
from source.util.ReflectionHandler import ReflectionHandler


class MatchedReferenceSequenceEncoder(DatasetEncoder):

    dataset_mapping = {
        "RepertoireDataset": "MatchedReferenceSequenceRepertoireEncoder"
    }

    def __init__(self, reference_sequences: list, summary: SequenceMatchingSummaryType, same_length_sequence: bool, metadata_attrs_to_match: list, max_edit_distance: int, chunk_size: int, name: str = None):

        self.reference_sequences = reference_sequences
        self.summary = summary
        self.same_length_sequence = same_length_sequence
        self.metadata_attrs_to_match = metadata_attrs_to_match
        self.max_edit_distance = max_edit_distance
        self.chunk_size = chunk_size
        self.name = name

    @staticmethod
    def _prepare_parameters(max_edit_distance: int, summary: str, reference_sequences: dict, metadata_attrs_to_match: list, same_length_sequence: bool, chunk_size: int = 128, name: str = None):

        location = "MatchedReferenceSequenceEncoder"

        ParameterValidator.assert_type_and_value(max_edit_distance, int, location, "max_edit_distance", min_inclusive=0)
        ParameterValidator.assert_keys(list(reference_sequences.keys()), ["format", "path", "params"], location, "reference_sequences", exclusive=False)
        ParameterValidator.assert_in_valid_list(summary.upper(), [item.name for item in SequenceMatchingSummaryType], location, "summary")

        valid_formats = ReflectionHandler.discover_classes_by_partial_name("SequenceImport", "sequence_import/")
        ParameterValidator.assert_in_valid_list(f"{reference_sequences['format']}SequenceImport", valid_formats, location,
                                                "format in reference_sequences")

        assert os.path.isfile(reference_sequences["path"]), f"{location}: the file {reference_sequences['path']} does not exist. " \
                                                            f"Specify the correct path under reference_sequences."

        importer = ReflectionHandler.get_class_by_name("{}SequenceImport".format(reference_sequences["format"]))

        sequences = importer.import_items(reference_sequences["path"], **reference_sequences.get("params", {}))

        return {
            "max_edit_distance": max_edit_distance,
            "summary": SequenceMatchingSummaryType[summary.upper()],
            "reference_sequences": sequences,
            "metadata_attrs_to_match": metadata_attrs_to_match,
            "same_length_sequence": same_length_sequence,
            "chunk_size": chunk_size,
            "name": name
        }

    @staticmethod
    def build_object(dataset=None, **params):
        try:
            prepared_parameters = MatchedReferenceSequenceEncoder._prepare_parameters(**params)
            encoder = ReflectionHandler.get_class_by_name(MatchedReferenceSequenceEncoder.dataset_mapping[dataset.__class__.__name__],
                                                          "reference_encoding/")(**prepared_parameters)
        except ValueError:
            raise ValueError("{} is not defined for dataset of type {}.".format(MatchedReferenceSequenceEncoder.__name__,
                                                                                dataset.__class__.__name__))
        return encoder

    def encode(self, dataset, params: EncoderParams):
        cache_key = CacheHandler.generate_cache_key(self._prepare_caching_params(dataset, params))
        print(cache_key)
        encoded_dataset = CacheHandler.memo(cache_key,
                                            lambda: self._encode_new_dataset(dataset, params))

        return encoded_dataset

    def _prepare_caching_params(self, dataset, params: EncoderParams):

        encoding_params_desc = {"same_length_sequence": self.same_length_sequence,
                                "summary": self.summary,
                                "metadata_attrs_to_match": sorted(self.metadata_attrs_to_match),
                                "max_edit_distance": self.max_edit_distance,
                                "reference_sequences": sorted([seq.get_sequence() for seq in self.reference_sequences])}

        return (("dataset_identifiers", tuple(dataset.get_example_ids())),
                ("dataset_metadata", dataset.metadata_file),
                ("dataset_type", dataset.__class__.__name__),
                ("labels", tuple(params["label_configuration"].get_labels_by_name())),
                ("encoding", MatchedReferenceSequenceEncoder.__name__),
                ("learn_model", params["learn_model"]),
                ("encoding_params", encoding_params_desc),)

    @abc.abstractmethod
    def _encode_new_dataset(self, dataset, params: EncoderParams):
        pass

    def store(self, encoded_dataset, params: EncoderParams):
        PickleExporter.export(encoded_dataset, params["result_path"])
