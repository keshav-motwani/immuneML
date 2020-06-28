from source.analysis.sequence_matching.SequenceMatchedRepertoire import SequenceMatchedRepertoire


class SequenceMatchedReferenceRepertoire(SequenceMatchedRepertoire):

    __slots__ = ['identifier', 'index', 'filename', 'metadata', 'chains', 'total_reads', 'unique_reads', 'reference_sequences']

    def __init__(self, identifier: str, index: int, filename: str, metadata: dict, chains: list, total_reads: int,
                 unique_reads: int, reference_sequences=None):
        super().__init__(identifier, index, filename, metadata, chains, total_reads, unique_reads)
        self.reference_sequences = reference_sequences
