class SequenceMatchedDataset(object):

    __slots__ = ['repertoires', 'reference_sequences']

    def __init__(self, repertoires: list, reference_sequences: list):
        self.repertoires = repertoires
        self.reference_sequences = reference_sequences
