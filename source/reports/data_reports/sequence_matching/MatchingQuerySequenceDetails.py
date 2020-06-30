import csv

from source.analysis.sequence_matching.query_sequence.QuerySequenceMatcher import QuerySequenceMatcher
from source.analysis.sequence_matching.query_sequence.SequenceMatchedQueryRepertoire import SequenceMatchedQueryRepertoire
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.encodings.reference_encoding.SequenceMatchingSummaryType import SequenceMatchingSummaryType
from source.reports.data_reports.DataReport import DataReport
from source.util.ParameterValidator import ParameterValidator
from source.util.PathBuilder import PathBuilder
from source.util.ReflectionHandler import ReflectionHandler


class MatchingQuerySequenceDetails(DataReport):
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

        location = "MatchingQuerySequenceDetails"

        if "max_edit_distance" in kwargs:
            ParameterValidator.assert_type_and_value(kwargs["max_edit_distance"], int, location, "max_edit_distance")

        if "reference_sequences" in kwargs:
            ParameterValidator.assert_keys(list(kwargs["reference_sequences"].keys()), ["format", "path", "params"], location, "reference_sequences", exclusive=False)

            importer = ReflectionHandler.get_class_by_name("{}SequenceImport".format(kwargs["reference_sequences"]["format"]))
            kwargs["reference_sequences"] = importer.import_items(kwargs["reference_sequences"]["path"], **kwargs["reference_sequences"].get("params", {})) \
                if kwargs["reference_sequences"] is not None else None

        return MatchingQuerySequenceDetails(**kwargs)


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
        self._make_overview()
        self._make_matching_report()

    def _match_repertoires(self):

        matcher = QuerySequenceMatcher()

        matched_info = matcher.match(dataset=self.dataset,
                                     reference_sequences=self.reference_sequences,
                                     same_length_sequence=self.same_length_sequence,
                                     metadata_attrs_to_match=self.metadata_attrs_to_match,
                                     max_edit_distance=self.max_edit_distance,
                                     batch_size=self.batch_size)

        self.matching_results = matched_info

        return matched_info

    def _make_overview(self):

        filename = self.result_path + "matching_sequence_overview.tsv"

        fieldnames = ["repertoire_identifier", "repertoire_size", "max_edit_distance", "metadata_attrs_to_match",
                      "same_length_sequence"] + [type.name.lower() for type in SequenceMatchingSummaryType]
        for label in self.dataset.params.keys():
            fieldnames.append("{}".format(label))

        self._write_rows(filename, fieldnames)

        return filename

    def _write_rows(self, filename: str, fieldnames: list):

        with open(filename, "w") as file:
            csv_writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter="\t")
            csv_writer.writeheader()

            for repertoire in self.matching_results.repertoires:

                row = {
                    "repertoire_identifier": repertoire.identifier,
                    "repertoire_size": len(repertoire.query_sequences),
                    "max_edit_distance": self.max_edit_distance,
                    "metadata_attrs_to_match": self.metadata_attrs_to_match,
                    "same_length_sequence": self.same_length_sequence
                }

                row = {
                    **row,
                    **{
                        name: getattr(repertoire, name) for name in [type.name.lower() for type in SequenceMatchingSummaryType]
                    }
                }

                for label in self.dataset.params.keys():
                    row[label] = repertoire.metadata[label]

                csv_writer.writerow(row)

    def _make_matching_report(self):

        filenames = []

        for repertoire in self.matching_results.repertoires:
            filenames.append(self._make_repertoire_report(repertoire))

        return filenames

    def _make_repertoire_report(self, repertoire: SequenceMatchedQueryRepertoire):

        filename = self.result_path + "{}_{}.tsv".format(repertoire.identifier,
                                                         list({repertoire.query_sequences[i].query_sequence["chain"]
                                                                        for i in range(len(repertoire.query_sequences))}))

        fieldnames = ["sequence", "v_gene", "j_gene", "chain", "clone_count", "matching_reference_sequences",
                      "max_edit_distance", "metadata_attrs_to_match", "same_length_sequence"]
        for label in self.dataset.params.keys():
            fieldnames.append("{}".format(label))

        with open(filename, "w") as file:

            csv_writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter="\t")
            csv_writer.writeheader()

            for sequence in repertoire.query_sequences:

                row = {
                    "sequence": sequence.query_sequence,
                    "v_gene": sequence.query_sequence["v_gene"],
                    "j_gene": sequence.query_sequence["j_gene"],
                    "chain": sequence.query_sequence["chain"],
                    "clone_count": sequence.query_sequence["count"],
                    "matching_reference_sequences": str(sequence.matching_reference_sequences),
                    "max_edit_distance": self.max_edit_distance,
                    "metadata_attrs_to_match": self.metadata_attrs_to_match,
                    "same_length_sequence": self.same_length_sequence
                }

                for label in self.dataset.params.keys():
                    row[label] = repertoire.metadata[label]

                csv_writer.writerow(row)

        return filename
