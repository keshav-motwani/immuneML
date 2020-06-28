from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence


class MatchedQuerySequence(object):

    __slots__ = ['query_sequence', 'matching_reference_sequences']

    def __init__(self, query_sequence, matching_reference_sequences):
        self.query_sequence = self.create_sequence_dict(query_sequence)
        self.matching_reference_sequences = [self.create_sequence_dict(sequence) for sequence in matching_reference_sequences]

    def create_sequence_dict(self, sequence: ReceptorSequence):

        return {"sequence": sequence.get_sequence(), **vars(sequence.metadata)}
