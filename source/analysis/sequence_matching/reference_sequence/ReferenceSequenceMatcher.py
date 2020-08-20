import time
from collections import defaultdict
from functools import partial
from multiprocessing.pool import Pool
import numpy as np

from source.analysis.sequence_matching.HashedReceptorSequence import HashedReceptorSequence
from source.analysis.sequence_matching.SequenceMatchedDataset import SequenceMatchedDataset
from source.analysis.sequence_matching.SequenceMatcher import SequenceMatcher
from source.analysis.sequence_matching.reference_sequence.MatchedReferenceSequence import MatchedReferenceSequence
from source.analysis.sequence_matching.reference_sequence.SequenceMatchedReferenceRepertoire import \
    SequenceMatchedReferenceRepertoire
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.repertoire.Repertoire import Repertoire


class ReferenceSequenceMatcher:

    @staticmethod
    def match(dataset: RepertoireDataset, reference_sequences: list, same_length_sequence: bool,
              metadata_attrs_to_match: list, max_edit_distance: int, batch_size: int) -> SequenceMatchedDataset:

        a = time.time()

        hashed_reference_list = [HashedReceptorSequence.hash_sequence(sequence, metadata_attrs_to_match) for
                                 sequence in reference_sequences]

        with Pool(batch_size, maxtasksperchild=1) as pool:
            fn = partial(ReferenceSequenceMatcher.match_repertoire, hashed_reference_list,
                         metadata_attrs_to_match, same_length_sequence, max_edit_distance)
            matched = pool.starmap(fn, enumerate(dataset.repertoires))

        b = time.time()

        print("Time elapsed in matching repertoires:", str(b - a))

        return SequenceMatchedDataset(matched, reference_sequences)

    @staticmethod
    def match_repertoire(hashed_reference_list: list,
                         metadata_attrs_to_match: list, same_length_sequence: bool,
                         max_edit_distance: int, index: int = 0,
                         repertoire: Repertoire = None) -> SequenceMatchedReferenceRepertoire:

        a = time.time()

        hashed_query_list = [HashedReceptorSequence.hash_sequence(sequence, metadata_attrs_to_match)
                             for sequence in repertoire.sequences if sequence.metadata.frame_type.upper() == "IN"]

        matches = SequenceMatcher.evaluate_repertoire_matches(hashed_query_list, hashed_reference_list,
                                                              same_length_sequence, max_edit_distance)

        hashes = np.array([hashed_query.hash for hashed_query in hashed_query_list])

        matches = ReferenceSequenceMatcher.generate_reference_to_query_map(matches, hashes)

        counts = repertoire.get_counts()
        total_reads = int(np.sum(counts))
        unique_reads = counts.shape[0]

        matched = SequenceMatchedReferenceRepertoire(
            identifier=repertoire.identifier,
            index=index,
            filename=repertoire.data_filename,
            metadata=repertoire.metadata,
            chains=None,
            total_reads=total_reads,
            unique_reads=unique_reads
        )

        matched.reference_sequences = [ReferenceSequenceMatcher.match_sequence(reference, hashed_query_list, matches)
                                       for reference in hashed_reference_list]

        b = time.time()

        print(repertoire.data_filename, "took", b - a, "seconds")

        return matched

    @staticmethod
    def match_sequence(hashed_reference: HashedReceptorSequence, hashed_query_list,
                       reference_to_query_matches_map: dict) -> MatchedReferenceSequence:

        matching_query_sequences = [hashed_query_list[i].sequence for i in
                                    reference_to_query_matches_map.get(hashed_reference.hash, [])]

        result = MatchedReferenceSequence(reference_sequence=hashed_reference.sequence,
                                          matching_query_sequences=matching_query_sequences)

        return result

    @staticmethod
    def generate_reference_to_query_map(matches: set, query_hashes):

        mapping = defaultdict(list)

        for pair in matches:
            mapping[pair[1]].append(np.nonzero(query_hashes == pair[0])[0])

        mapping = {i: np.hstack(j) for i, j in mapping.items()}

        return mapping