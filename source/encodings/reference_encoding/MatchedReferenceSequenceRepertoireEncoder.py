import math
from multiprocessing.pool import Pool

import numpy as np
import pandas as pd
from scipy import sparse

from source.analysis.sequence_matching.HashedReceptorSequence import HashedReceptorSequence
from source.analysis.sequence_matching.reference_sequence.ReferenceSequenceMatcher import ReferenceSequenceMatcher
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.data_model.repertoire.Repertoire import Repertoire
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.MatchedReferenceSequenceEncoder import MatchedReferenceSequenceEncoder
from source.encodings.reference_encoding.SequenceMatchingSummaryType import SequenceMatchingSummaryType


class MatchedReferenceSequenceRepertoireEncoder(MatchedReferenceSequenceEncoder):

    def _encode_new_dataset(self, dataset, params: EncoderParams):

        self.hashed_reference_sequences = self._hash_reference_sequences()

        encoded_repertoires, repertoire_names, feature_names, feature_annotations, labels = \
            self._encode_repertoires(dataset, params)

        encoded_dataset = RepertoireDataset(repertoires=dataset.repertoires, params=dataset.params,
                                            metadata_file=dataset.metadata_file)

        encoded_dataset.add_encoded_data(EncodedData(
            examples=encoded_repertoires,
            labels=labels,
            feature_names=feature_names,
            feature_annotations=feature_annotations,
            example_ids=repertoire_names,
            encoding=MatchedReferenceSequenceEncoder.__name__
        ))

        self.store(encoded_dataset, params)
        return encoded_dataset

    def _encode_repertoires(self, dataset, params: EncoderParams):

        features = [self.get_feature_name(sequence, self.metadata_attrs_to_match) for sequence in self.reference_sequences]

        feature_annotations = []

        for sequence in self.reference_sequences:
            feature_annotation = {
                **{k: v for k, v in vars(sequence.metadata).items() if k != "custom_params"},
                **sequence.metadata.custom_params
            }
            feature_annotations.append(feature_annotation)

        feature_annotations = pd.DataFrame(feature_annotations)

        arguments = [(params, index, repertoire) for index, repertoire in enumerate(dataset.repertoires)]

        with Pool(params["batch_size"]) as pool:
            chunksize = math.floor(dataset.get_example_count()/params["batch_size"]) + 1
            repertoires = pool.starmap(self._encode_repertoire, arguments, chunksize=chunksize)

        encoded_repertoire_list, repertoire_names, labels = zip(*repertoires)

        encoded_labels = {k: [dic[k] for dic in labels] for k in labels[0]}

        encoded_repertoires = sparse.csr_matrix(np.vstack(encoded_repertoire_list))

        return encoded_repertoires, repertoire_names, features, feature_annotations, encoded_labels

    def _encode_repertoire(self, params: EncoderParams, index, repertoire: Repertoire):

        encoded = np.zeros((len(self.reference_sequences,)),
                                                 dtype=float)

        rep = ReferenceSequenceMatcher.match_repertoire(self.hashed_reference_sequences,
                                                  self.metadata_attrs_to_match,
                                                  self.same_length_sequence,
                                                  self.max_edit_distance,
                                                  index,
                                                  repertoire)

        for j in range(len(self.reference_sequences)):
            encoded[j] = self.get_feature_value(rep, j)

        label_config = params["label_configuration"]
        labels = dict()

        for label_name in label_config.get_labels_by_name():
            label = repertoire.metadata[label_name]
            labels[label_name] = label

        repertoire.free_memory()

        return encoded, repertoire.identifier, labels

    def get_feature_name(self, sequence, metadata_attrs):
        if len(metadata_attrs) == 0:
            feature = sequence.get_sequence()
        else:
            feature = sequence.get_sequence() + "-" + "-".join([getattr(sequence.metadata, k) for k in metadata_attrs])
        return feature

    def get_feature_value(self, matched_repertoire, j):

        if self.summary == SequenceMatchingSummaryType.UNIQUE_READS_WITH_MATCH:
            value = len(matched_repertoire.reference_sequences[j].matching_query_sequences)

        elif self.summary == SequenceMatchingSummaryType.PCT_UNIQUE_READS_WITH_MATCH:
            value = 100 * len(matched_repertoire.reference_sequences[j].matching_query_sequences) / \
                    matched_repertoire.unique_reads

        elif self.summary == SequenceMatchingSummaryType.TOTAL_READS_WITH_MATCH:
            value = sum([metadata["count"] for metadata in
                         matched_repertoire.reference_sequences[j].matching_query_sequences])

        else:
            value = 100 * sum([metadata["count"] for metadata in
                         matched_repertoire.reference_sequences[j].matching_query_sequences]) / \
                    matched_repertoire.total_reads

        return value

    def _hash_reference_sequences(self):

        hashed_reference_list = [HashedReceptorSequence.hash_sequence(sequence, self.metadata_attrs_to_match) for
                                 sequence in self.reference_sequences]

        return hashed_reference_list
