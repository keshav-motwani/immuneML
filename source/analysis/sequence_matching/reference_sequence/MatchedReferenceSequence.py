from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence


class MatchedReferenceSequence(object):

    __slots__ = ['reference_sequence', 'matching_query_sequences']

    def __init__(self, reference_sequence, matching_query_sequences):
        self.reference_sequence = self.create_sequence_dict(reference_sequence)
        self.matching_query_sequences = [self.create_sequence_dict(sequence) for sequence in matching_query_sequences]

    def create_sequence_dict(self, sequence: ReceptorSequence):

        return {"sequence": sequence.get_sequence(), **vars(sequence.metadata), "identifier": sequence.identifier}
