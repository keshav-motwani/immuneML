import time
from collections import defaultdict
from functools import partial
from multiprocessing.pool import Pool

from source.analysis.sequence_matching.HashedReceptorSequence import HashedReceptorSequence
from source.analysis.sequence_matching.SequenceMatchedDataset import SequenceMatchedDataset
from source.analysis.sequence_matching.SequenceMatcher import SequenceMatcher
from source.analysis.sequence_matching.query_sequence.MatchedQuerySequence import MatchedQuerySequence
from source.analysis.sequence_matching.query_sequence.SequenceMatchedQueryRepertoire import SequenceMatchedQueryRepertoire
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.repertoire.Repertoire import Repertoire


class QuerySequenceMatcher:

    @staticmethod
    def match(dataset: RepertoireDataset, reference_sequences: list, same_length_sequence: bool,
              metadata_attrs_to_match: list, max_edit_distance: int, batch_size: int,
              return_sequence_info=True) -> SequenceMatchedDataset:

        a = time.time()

        hashed_reference_list = [HashedReceptorSequence.hash_sequence(sequence, metadata_attrs_to_match) for
                                 sequence in reference_sequences]

        with Pool(batch_size, maxtasksperchild=1) as pool:
            fn = partial(QuerySequenceMatcher.match_repertoire, hashed_reference_list,
                         metadata_attrs_to_match, same_length_sequence, max_edit_distance,
                         return_sequence_info)
            matched = pool.starmap(fn, enumerate(dataset.repertoires), chunksize=1)

        b = time.time()

        print("Time elapsed in matching repertoires:", str(b - a))

        return SequenceMatchedDataset(matched, reference_sequences)

    @staticmethod
    def match_repertoire(hashed_reference_list: list,
                         metadata_attrs_to_match: list, same_length_sequence: bool,
                         max_edit_distance: int, return_sequence_info=True,
                         index: int = 0, repertoire: Repertoire = None) -> SequenceMatchedQueryRepertoire:

        a = time.time()

        hashed_query_list = [HashedReceptorSequence.hash_sequence(sequence, metadata_attrs_to_match)
                             for sequence in repertoire.sequences if sequence.metadata.frame_type.upper() == "IN"]

        matches = SequenceMatcher.evaluate_repertoire_matches(hashed_query_list, hashed_reference_list,
                                                              same_length_sequence, max_edit_distance)

        matches = QuerySequenceMatcher.generate_query_to_reference_map(matches)

        total_reads = sum([sequence.metadata.count for sequence in repertoire.sequences])
        unique_reads = len(repertoire.sequences)

        matched = SequenceMatchedQueryRepertoire(
            identifier=repertoire.identifier,
            index=index,
            filename=repertoire.data_filename,
            metadata=repertoire.metadata,
            chains=list(set([sequence.metadata.chain for sequence in repertoire.sequences])),
            total_reads=total_reads,
            unique_reads=unique_reads
        )

        matched_query_sequences = [QuerySequenceMatcher.match_sequence(hashed_query, hashed_reference_list, matches)
                                   for hashed_query in hashed_query_list]

        matched.total_reads_with_match = QuerySequenceMatcher.compute_total_reads_with_match(matched_query_sequences)
        matched.pct_total_reads_with_match = matched.total_reads_with_match / total_reads
        matched.unique_reads_with_match = QuerySequenceMatcher.compute_unique_reads_with_match(matched_query_sequences)
        matched.pct_unique_reads_with_match = matched.unique_reads_with_match / unique_reads

        if return_sequence_info:
            matched.query_sequences = matched_query_sequences

        b = time.time()

        print(repertoire.data_filename, "took", b - a, "seconds")

        return matched

    @staticmethod
    def match_sequence(hashed_query: HashedReceptorSequence, hashed_reference_list: list,
                       query_to_reference_matches_map: dict) -> MatchedQuerySequence:

        if hashed_query.hash in query_to_reference_matches_map:
            matching_reference_sequences = [hashed_reference.sequence for hashed_reference in hashed_reference_list if
                                            hashed_reference.hash in query_to_reference_matches_map[hashed_query.hash]]
        else:
            matching_reference_sequences = []

        result = MatchedQuerySequence(query_sequence=hashed_query.sequence,
                                      matching_reference_sequences=matching_reference_sequences)

        return result

    @staticmethod
    def compute_total_reads_with_match(matched_query_sequences):

        total_reads = sum([matched_query.query_sequence["count"] for matched_query in matched_query_sequences
                           if len(matched_query.matching_reference_sequences) > 0])

        return total_reads

    @staticmethod
    def compute_unique_reads_with_match(matched_query_sequences):

        unique_reads = len([matched_query.query_sequence["count"] for matched_query in matched_query_sequences
                            if len(matched_query.matching_reference_sequences) > 0])

        return unique_reads

    @staticmethod
    def generate_query_to_reference_map(matches: set):

        mapping = defaultdict(set)

        for pair in matches:
            mapping[pair[0]].add(pair[1])

        return mapping
