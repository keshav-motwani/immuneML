import os.path

import pandas as pd
from source.IO.dataset_export.PickleExporter import PickleExporter
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.encodings.DatasetEncoder import DatasetEncoder
from source.encodings.EncoderParams import EncoderParams
from source.util.PathBuilder import PathBuilder


class DeepRCEncoder(DatasetEncoder):
    """
    DeepRCEncoder should be used in combination with the DeepRC ML method (:py:mod:`source.ml_methods.DeepRC.DeepRC`).
    This encoder writes the data in a RepertoireDataset to .tsv files.
    For each repertoire, one .tsv file is created containing the amino acid sequences and the counts.
    Additionally one metadata .tsv file is created, which describes the subset of repertoires that is encoded by
    a given instance of the DeepRCEncoder.

    Specification:

    .. indent with spaces
    .. code-block:: yaml

        my_deeprc_encoder: DeepRC
    """
    ID_COLUMN = "ID"
    SEQUENCE_COLUMN = "amino_acid"
    COUNTS_COLUMN = "templates"
    SEP = "\t"
    EXTENSION = "tsv"

    def __init__(self, context: dict = None, name: str = None):
        self.context = context
        self.name = name

    def set_context(self, context: dict):
        self.context = context
        return self

    @staticmethod
    def _prepare_parameters(context: dict = None, name: str = None):
        return {"context": context,
                "name": name}

    @staticmethod
    def build_object(dataset, **params):
        if isinstance(dataset, RepertoireDataset):
            prepared_params = DeepRCEncoder._prepare_parameters(**params)
            return DeepRCEncoder(**prepared_params)
        else:
            raise ValueError("DeepRCEncoder is not defined for dataset types which are not RepertoireDataset.")

    def export_repertoire_tsv_files(self, output_folder):
        repertoires = self.context["dataset"].repertoires

        for repertoire in repertoires:
            filepath = f"{output_folder}/{repertoire.identifier}.{DeepRCEncoder.EXTENSION}"

            if not os.path.isfile(filepath):
                df = pd.DataFrame({DeepRCEncoder.SEQUENCE_COLUMN: repertoire.get_sequence_aas(), DeepRCEncoder.COUNTS_COLUMN: repertoire.get_counts()})
                df.to_csv(path_or_buf=filepath, sep=DeepRCEncoder.SEP, index=False)

    def export_metadata_file(self, dataset, labels, output_folder):
        metadata_filepath = f"{output_folder}/{dataset.identifier}_metadata.{DeepRCEncoder.EXTENSION}"
        metadata = dataset.get_metadata(labels, return_df=True)
        metadata[DeepRCEncoder.ID_COLUMN] = dataset.get_repertoire_ids()

        metadata.to_csv(path_or_buf=metadata_filepath, sep="\t", index=False)

        return metadata_filepath

    def encode(self, dataset, params: EncoderParams) -> RepertoireDataset:
        result_path = params["result_path"] + "/encoding"
        PathBuilder.build(result_path)

        self.export_repertoire_tsv_files(result_path)

        labels = params["label_configuration"].get_labels_by_name()
        metadata_filepath = self.export_metadata_file(dataset, labels, result_path)

        encoded_dataset = dataset.clone()
        encoded_dataset.encoded_data = EncodedData(examples=None, labels=dataset.get_metadata(labels),
                                                   example_ids=dataset.repertoire_ids,
                                                   encoding=DeepRCEncoder.__name__,
                                                   info={"metadata_filepath": metadata_filepath})

        self.store(encoded_dataset, params)

        return encoded_dataset


    def store(self, encoded_dataset, params: EncoderParams):
        PickleExporter.export(encoded_dataset, params["result_path"])