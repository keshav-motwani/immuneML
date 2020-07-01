import numpy as np
import pandas as pd
from scipy import sparse

from source.analysis.sequence_matching.reference_sequence.ReferenceSequenceMatcher import ReferenceSequenceMatcher
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.MatchedReferenceSequenceEncoder import MatchedReferenceSequenceEncoder
from source.encodings.reference_encoding.SequenceMatchingSummaryType import SequenceMatchingSummaryType


class MatchedReferenceSequenceRepertoireEncoder(MatchedReferenceSequenceEncoder):

    def _encode_new_dataset(self, dataset, params: EncoderParams):

        matched_info = self._match_repertoires(dataset, params["batch_size"])

        encoded_dataset = RepertoireDataset(repertoires=dataset.repertoires, params=dataset.params,
                                            metadata_file=dataset.metadata_file)

        encoded_repertoires, feature_names, feature_annotations, labels = self._encode_repertoires(dataset,
                                                                                                   matched_info, params)

        encoded_dataset.add_encoded_data(EncodedData(
            examples=encoded_repertoires,
            labels=labels,
            feature_names=feature_names,
            feature_annotations=feature_annotations,
            example_ids=[repertoire.identifier for repertoire in matched_info.repertoires],
            encoding=MatchedReferenceSequenceEncoder.__name__
        ))

        self.store(encoded_dataset, params)
        return encoded_dataset

    def _encode_repertoires(self, dataset, matched_info, params: EncoderParams):

        encoded_repertoires = np.zeros((dataset.get_example_count(), len(matched_info.reference_sequences)),
                                       dtype=float)

        labels = {label: [] for label in params["label_configuration"].get_labels_by_name()}

        features = [self.get_feature_name(sequence, self.metadata_attrs_to_match) for sequence in matched_info.reference_sequences]

        feature_annotations = []

        for i in range(len(matched_info.repertoires)):

            for j in range(len(features)):
                assert matched_info.reference_sequences[j].get_sequence() == matched_info.repertoires[i].reference_sequences[
                    j].reference_sequence["sequence"], "Order of reference_sequences is not consistent across repertoires."
                encoded_repertoires[i, j] = self.get_feature_value(matched_info, i, j)

            for label_index, label in enumerate(params["label_configuration"].get_labels_by_name()):
                labels[label].append(matched_info.repertoires[i].metadata[label])

        for sequence in matched_info.reference_sequences:
            feature_annotation = {
                **{k: v for k, v in vars(sequence.metadata).items() if k != "custom_params"},
                **sequence.metadata.custom_params
            }
            feature_annotations.append(feature_annotation)

        feature_annotations = pd.DataFrame(feature_annotations)

        encoded_repertoires = sparse.csr_matrix(encoded_repertoires)

        return encoded_repertoires, [str(feature) for feature in features], feature_annotations, labels

    def get_feature_name(self, sequence, metadata_attrs):
        return (sequence.get_sequence() + str(tuple([getattr(sequence.metadata, k) for k in metadata_attrs]))).replace("\"\"", "")

    def get_feature_value(self, matched_info, i, j):

        if self.summary == SequenceMatchingSummaryType.UNIQUE_READS_WITH_MATCH:
            value = len(matched_info.repertoires[i].reference_sequences[j].matching_query_sequences)

        elif self.summary == SequenceMatchingSummaryType.PCT_UNIQUE_READS_WITH_MATCH:
            value = len(matched_info.repertoires[i].reference_sequences[j].matching_query_sequences) / \
                    matched_info.repertoires[i].unique_reads

        elif self.summary == SequenceMatchingSummaryType.TOTAL_READS_WITH_MATCH:
            value = sum([metadata["count"] for metadata in
                         matched_info.repertoires[i].reference_sequences[j].matching_query_sequences])

        else:
            value = sum([metadata["count"] for metadata in
                         matched_info.repertoires[i].reference_sequences[j].matching_query_sequences]) / \
                    matched_info.repertoires[i].total_reads

        return value

    def _match_repertoires(self, dataset: RepertoireDataset, batch_size: int):

        matcher = ReferenceSequenceMatcher

        matched_info = matcher.match(dataset=dataset,
                                     reference_sequences=self.reference_sequences,
                                     max_edit_distance=self.max_edit_distance,
                                     same_length_sequence=self.same_length_sequence,
                                     metadata_attrs_to_match=self.metadata_attrs_to_match,
                                     batch_size=batch_size)

        return matched_info
