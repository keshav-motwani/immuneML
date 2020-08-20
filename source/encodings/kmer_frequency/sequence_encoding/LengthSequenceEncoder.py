import warnings

from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.data_model.receptor.receptor_sequence.SequenceFrameType import SequenceFrameType
from source.encodings.EncoderParams import EncoderParams
from source.encodings.kmer_frequency.sequence_encoding.SequenceEncodingStrategy import SequenceEncodingStrategy
from source.util.KmerHelper import KmerHelper


class LengthSequenceEncoder(SequenceEncodingStrategy):

    @staticmethod
    def encode_sequence(sequence: ReceptorSequence, params: EncoderParams):
        if sequence.metadata is not None and sequence.metadata.frame_type.upper() != SequenceFrameType.IN.name:
            warnings.warn('Sequence either has out or stop codon. Ignoring sequence.')
            return None
        else:
            length = len(sequence.get_sequence())
            return [str(length)]

    @staticmethod
    def get_feature_names(params: EncoderParams):
        return ["length"]
