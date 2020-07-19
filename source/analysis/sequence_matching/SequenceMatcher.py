from collections import defaultdict

import Levenshtein_search
import distance

from source.analysis.sequence_matching.SlimReceptorSequence import SlimReceptorSequence
from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence


class SequenceMatcher:

    distance_function = None

    @staticmethod
    def generate_slim_sequence_map(hashed_sequence_list: list):

        mapping = defaultdict(set)

        for sequence in hashed_sequence_list:
            mapping[sequence.slim_sequence.sequence].add(sequence.slim_sequence)

        return mapping

    @staticmethod
    def evaluate_repertoire_matches(hashed_query_list: list, hashed_reference_list: list, same_length_sequence: bool, max_edit_distance: int):

        slim_query_sequence_map = SequenceMatcher.generate_slim_sequence_map(hashed_query_list)

        repertoire_sequences = list(slim_query_sequence_map.keys())
        wordset = Levenshtein_search.populate_wordset(-1, repertoire_sequences)

        result = []

        for i in hashed_reference_list:

            string_query_matches = Levenshtein_search.lookup(wordset, i.slim_sequence.sequence, max_edit_distance)
            string_query_matches = [sequence[0] for sequence in string_query_matches]

            slim_query_matches = [slim_query_sequence_map[sequence] for sequence in string_query_matches]
            slim_query_matches = [item for items in slim_query_matches for item in items if SequenceMatcher.evaluate_slim_sequence_match(i.slim_sequence, item, same_length_sequence)]

            tmp = [(hash(j), i.hash) for j in slim_query_matches]

            result.extend(tmp)

        return result

    @staticmethod
    def evaluate_slim_sequence_match(slim_sequence_1: SlimReceptorSequence, slim_sequence_2: SlimReceptorSequence,
                                     same_length_sequence: bool):

        match = True

        if same_length_sequence:
            match = len(slim_sequence_1.sequence) == len(slim_sequence_2.sequence)

        if match:
            match = slim_sequence_1.metadata == slim_sequence_2.metadata

        return match

    @staticmethod
    def matches_sequence(original_sequence: ReceptorSequence, reference_sequence: ReceptorSequence, max_distance):
        """
        :param original_sequence: ReceptorSequence
        :param reference_sequence: ReceptorSequence
        :param max_distance: max allowed Levenshtein distance between two sequences to be considered a match
        :return: True if chain, v_gene and j_gene are the same and sequences are within given Levenshtein distance
        """
        return reference_sequence.metadata.chain == original_sequence.metadata.chain \
            and SequenceMatcher.matches_gene(reference_sequence.metadata.v_gene, original_sequence.metadata.v_gene) \
            and SequenceMatcher.matches_gene(reference_sequence.metadata.j_gene, original_sequence.metadata.j_gene) \
            and distance.levenshtein(original_sequence.get_sequence(), reference_sequence.get_sequence()) <= max_distance

    @staticmethod
    def matches_gene(gene1, gene2):
        if gene1 == gene2:
            return True
        else:
            return gene2.split("-", 1)[0] == gene1 or gene1.split("-", 1)[0] == gene2
