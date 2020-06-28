from enum import Enum


class SequenceMatchingSummaryType(Enum):

    COUNT = 0
    CLONAL_PERCENTAGE = 1
    PERCENTAGE = 2

    # total number of reads of query sequences that match some reference sequence
    # use case: interested in number of reads in a repertoire that map to a set of reference sequences
    TOTAL_READS_WITH_MATCH = 3

    # total number of reads of query sequences that match some reference sequence normalized by the total number of
    #   reads in the query sequences
    # use case: interested in percentage of reads in a repertoire that map to a set of reference sequences
    PCT_TOTAL_READS_WITH_MATCH = 4

    # number of query sequences that match some reference sequence
    # use case: interested in number of unique sequences in a repertoire that map to a set of reference sequences
    UNIQUE_READS_WITH_MATCH = 5

    # number of query sequences that match some reference sequence normalized by the number of query sequences
    # use case: interested in percentage of unique sequences in a repertoire that map to a set of reference sequences
    PCT_UNIQUE_READS_WITH_MATCH = 6
