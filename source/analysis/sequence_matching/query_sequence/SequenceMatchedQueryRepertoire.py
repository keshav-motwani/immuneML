from source.analysis.sequence_matching.SequenceMatchedRepertoire import SequenceMatchedRepertoire
from source.encodings.reference_encoding.SequenceMatchingSummaryType import SequenceMatchingSummaryType


class SequenceMatchedQueryRepertoire(SequenceMatchedRepertoire):
    __slots__ = ['query_sequences'] + [summary_type.name.lower() for summary_type in SequenceMatchingSummaryType]

    def __init__(self, identifier: str, index: int, filename: str, metadata: dict, chains: list, total_reads: int,
                 unique_reads: int, query_sequences=None):
        super().__init__(identifier, index, filename, metadata, chains, total_reads, unique_reads)
        self.query_sequences = query_sequences

        for summary_type in SequenceMatchingSummaryType:
            setattr(self, summary_type.name.lower(), None)
