import csv

from source.analysis.sequence_matching.reference_sequence.ReferenceSequenceMatcher import ReferenceSequenceMatcher
from source.analysis.sequence_matching.reference_sequence.SequenceMatchedReferenceRepertoire import \
    SequenceMatchedReferenceRepertoire
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.reports.data_reports.DataReport import DataReport
from source.util.ParameterValidator import ParameterValidator
from source.util.PathBuilder import PathBuilder
from source.util.ReflectionHandler import ReflectionHandler


class MatchingReferenceSequenceDetails(DataReport):
    """
    TODO: write description here
    params:
        - list of reference sequences
        - max Levenshtein distance
        - summary:  * count the number of sequences from the repertoire matched,
                    * get the percentage of sequences from the repertoire matched,
                    * get the percentage of sequences from the repertoire matched with respect to clonal counts
    """

    @classmethod
    def build_object(cls, **kwargs):

        location = "MatchingSequenceDetails"

        if "max_edit_distance" in kwargs:
            ParameterValidator.assert_type_and_value(kwargs["max_edit_distance"], int, location, "max_edit_distance")

        if "reference_sequences" in kwargs:
            ParameterValidator.assert_keys(list(kwargs["reference_sequences"].keys()), ["format", "path"], location, "reference_sequences")

            importer = ReflectionHandler.get_class_by_name("{}SequenceImport".format(kwargs["reference_sequences"]["format"]))
            kwargs["reference_sequences"] = importer.import_items(kwargs["reference_sequences"]["path"]) \
                if kwargs["reference_sequences"] is not None else None

        return MatchingReferenceSequenceDetails(**kwargs)

    def __init__(self, dataset: RepertoireDataset = None, max_edit_distance: int = None,
                 metadata_attrs_to_match: list = [], same_length_sequence: bool = False,
                 reference_sequences: list = None, batch_size: int = 1, result_path: str = None):
        DataReport.__init__(self, dataset=dataset, result_path=result_path)
        self.dataset = dataset
        self.max_edit_distance = max_edit_distance
        self.metadata_attrs_to_match = metadata_attrs_to_match
        self.same_length_sequence = same_length_sequence
        self.reference_sequences = reference_sequences
        self.matching_results = None
        self.result_path = result_path
        self.batch_size = batch_size

    def generate(self):
        PathBuilder.build(self.result_path)
        self._match_repertoires()
        self._make_matching_report()

    def check_prerequisites(self):
        pass

    def _match_repertoires(self):

        matcher = ReferenceSequenceMatcher()

        matched_info = matcher.match(dataset=self.dataset,
                                     reference_sequences=self.reference_sequences,
                                     same_length_sequence=self.same_length_sequence,
                                     metadata_attrs_to_match=self.metadata_attrs_to_match,
                                     max_edit_distance=self.max_edit_distance,
                                     batch_size=self.batch_size)

        self.matching_results = matched_info

        return matched_info

    def _make_matching_report(self):

        filenames = []

        for repertoire in self.matching_results:
            filenames.append(self._make_repertoire_report(repertoire))

        return filenames

    def _make_repertoire_report(self, repertoire: SequenceMatchedReferenceRepertoire):

        filename = self.result_path + "{}.tsv".format(repertoire.identifier)

        with open(filename, "w") as file:
            csv_writer = csv.DictWriter(file,
                                        fieldnames=["sequence", "v_gene", "j_gene", "chain", "matching_query_sequences",
                                                    "max_edit_distance", "metadata_attrs_to_match", "same_length_sequence"],
                                        delimiter="\t")
            csv_writer.writeheader()

            for sequence in repertoire.reference_sequences:
                csv_writer.writerow({
                    "sequence": sequence.reference_sequence,
                    "v_gene": sequence.reference_sequence["v_gene"],
                    "j_gene": sequence.reference_sequence["j_gene"],
                    "chain": sequence.reference_sequence["chain"],
                    "matching_query_sequences": str(sequence.matching_query_sequences),
                    "max_edit_distance": self.max_edit_distance,
                    "metadata_attrs_to_match": self.metadata_attrs_to_match,
                    "same_length_sequence": self.same_length_sequence
                })

        return filename
