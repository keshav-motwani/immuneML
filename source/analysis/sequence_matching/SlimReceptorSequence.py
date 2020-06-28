from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence


class SlimReceptorSequence(object):

    __slots__ = ['sequence', 'metadata']

    def __init__(self, sequence, metadata):
        self.sequence = sequence
        self.metadata = metadata

    def __eq__(self, other):
        return self.sequence == other.sequence and self.metadata == other.metadata

    def __hash__(self):
        return hash((self.sequence, self.metadata))

    def __str__(self):
        return self.sequence + str(self.metadata)

    @staticmethod
    def slim_sequence(sequence: ReceptorSequence, metadata_attrs: list):
        return SlimReceptorSequence(sequence.get_sequence(),
                                    hash(tuple([getattr(sequence.metadata, k) for k in metadata_attrs])))
