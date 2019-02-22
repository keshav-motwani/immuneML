import abc

from source.data_model.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.encodings.EncoderParams import EncoderParams


class SequenceEncodingStrategy(metaclass=abc.ABCMeta):

    @staticmethod
    @abc.abstractmethod
    def encode_sequence(sequence: ReceptorSequence, params: EncoderParams):
        pass
