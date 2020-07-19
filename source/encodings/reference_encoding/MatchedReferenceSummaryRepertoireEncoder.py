import math
from multiprocessing.pool import Pool

import numpy as np
from scipy import sparse
import pandas as pd

from source.analysis.sequence_matching.HashedReceptorSequence import HashedReceptorSequence
from source.analysis.sequence_matching.query_sequence.QuerySequenceMatcher import QuerySequenceMatcher
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.data_model.repertoire.Repertoire import Repertoire
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.MatchedReferenceSummaryEncoder import MatchedReferenceSummaryEncoder


class MatchedReferenceSummaryRepertoireEncoder(MatchedReferenceSummaryEncoder):

    def _encode_new_dataset(self, dataset, params: EncoderParams):

        self.hashed_reference_sequences = self._hash_reference_sequences()

        encoded_repertoires, repertoire_names, labels = self._encode_repertoires(dataset, params)

        encoded_dataset = RepertoireDataset(repertoires=dataset.repertoires, params=dataset.params,
                                            metadata_file=dataset.metadata_file)

        feature_names = [summary.name.lower() for summary in self.summary]

        encoded_dataset.add_encoded_data(EncodedData(
            examples=sparse.csr_matrix(encoded_repertoires),
            labels=labels,
            feature_names=feature_names,
            feature_annotations=pd.DataFrame({"feature": feature_names}),
            example_ids=repertoire_names,
            encoding=MatchedReferenceSummaryEncoder.__name__
        ))

        self.store(encoded_dataset, params)
        return encoded_dataset

    def _encode_repertoires(self, dataset, params: EncoderParams):

        arguments = [(params, index, repertoire) for index, repertoire in enumerate(dataset.repertoires)]

        with Pool(params["batch_size"]) as pool:
            chunksize = math.floor(dataset.get_example_count() / params["batch_size"]) + 1
            repertoires = pool.starmap(self._encode_repertoire, arguments, chunksize=chunksize)

        encoded_repertoire_list, repertoire_names, labels = zip(*repertoires)

        encoded_labels = {k: [dic[k] for dic in labels] for k in labels[0]}

        encoded_repertoires = sparse.csr_matrix(np.vstack(encoded_repertoire_list))

        return encoded_repertoires, repertoire_names, encoded_labels

    def _encode_repertoire(self, params: EncoderParams, index, repertoire: Repertoire):

        encoded = np.zeros((len(self.summary, )), dtype=float)

        rep = QuerySequenceMatcher.match_repertoire(self.hashed_reference_sequences,
                                                    self.metadata_attrs_to_match,
                                                    self.same_length_sequence,
                                                    self.max_edit_distance,
                                                    False,
                                                    index,
                                                    repertoire)

        for j, feature in enumerate(self.summary):
            encoded[j] = getattr(rep, feature.name.lower())

        label_config = params["label_configuration"]
        labels = dict()

        for label_name in label_config.get_labels_by_name():
            label = repertoire.metadata[label_name]
            labels[label_name] = label

        repertoire.free_memory()

        return encoded, repertoire.identifier, labels

    def _hash_reference_sequences(self):

        hashed_reference_list = [HashedReceptorSequence.hash_sequence(sequence, self.metadata_attrs_to_match) for
                                 sequence in self.reference_sequences]

        return hashed_reference_list
