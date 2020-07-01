import numpy as np
from scipy import sparse
import pandas as pd

from source.analysis.sequence_matching.query_sequence.QuerySequenceMatcher import QuerySequenceMatcher
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.MatchedReferenceSummaryEncoder import MatchedReferenceSummaryEncoder


class MatchedReferenceSummaryRepertoireEncoder(MatchedReferenceSummaryEncoder):

    def _encode_new_dataset(self, dataset, params: EncoderParams):

        matched_info = self._match_repertoires(dataset, params["batch_size"])

        encoded_dataset = RepertoireDataset(repertoires=dataset.repertoires,
                                            params=dataset.params,
                                            metadata_file=dataset.metadata_file)

        encoded_repertoires, labels = self._encode_repertoires(dataset, matched_info, params)

        feature_names = [summary.name.lower() for summary in self.summary]

        encoded_dataset.add_encoded_data(EncodedData(
            examples=sparse.csr_matrix(encoded_repertoires),
            labels=labels,
            feature_names=feature_names,
            feature_annotations=pd.DataFrame({"feature": feature_names}),
            example_ids=[repertoire.identifier for repertoire in matched_info.repertoires],
            encoding=MatchedReferenceSummaryEncoder.__name__
        ))

        self.store(encoded_dataset, params)
        return encoded_dataset

    def _encode_repertoires(self, dataset, matched_info, params: EncoderParams):

        encoded_repertoires = np.zeros((dataset.get_example_count(), len(self.summary)), dtype=float)

        labels = {label: [] for label in params["label_configuration"].get_labels_by_name()}

        for i in range(len(matched_info.repertoires)):

            for j, feature in enumerate(self.summary):
                encoded_repertoires[i, j] = getattr(matched_info.repertoires[i], feature.name.lower())

            for label_index, label in enumerate(params["label_configuration"].get_labels_by_name()):
                labels[label].append(matched_info.repertoires[i].metadata[label])

        return encoded_repertoires, labels

    def _match_repertoires(self, dataset: RepertoireDataset, batch_size: int):
        matcher = QuerySequenceMatcher
        matched_info = matcher.match(dataset=dataset,
                                     reference_sequences=self.reference_sequences,
                                     same_length_sequence=self.same_length_sequence,
                                     metadata_attrs_to_match=self.metadata_attrs_to_match,
                                     max_edit_distance=self.max_edit_distance,
                                     batch_size=batch_size,
                                     return_sequence_info=False)

        return matched_info
