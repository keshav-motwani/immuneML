class SequenceMatchedRepertoire(object):

    __slots__ = ['identifier', 'index', 'filename', 'metadata', 'chains', 'total_reads', 'unique_reads']

    def __init__(self, identifier: str, index: int, filename: str, metadata: dict, chains: list, total_reads: int,
                 unique_reads: int):
        self.identifier = identifier
        self.index = index
        self.filename = filename
        self.metadata = metadata
        self.chains = chains
        self.total_reads = total_reads
        self.unique_reads = unique_reads
