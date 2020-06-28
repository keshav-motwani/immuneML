from source.analysis.sequence_matching.SlimReceptorSequence import SlimReceptorSequence
from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence


class HashedReceptorSequence:

    __slots__ = ['hash', 'sequence', 'slim_sequence']

    def __init__(self, hash: int, sequence: ReceptorSequence, slim_sequence: SlimReceptorSequence):
        self.hash = hash
        self.sequence = sequence
        self.slim_sequence = slim_sequence

    @staticmethod
    def hash_sequence(sequence: ReceptorSequence, metadata_attrs):
        slim_sequence = SlimReceptorSequence.slim_sequence(sequence, metadata_attrs)
        return HashedReceptorSequence(hash(slim_sequence), sequence, slim_sequence)
