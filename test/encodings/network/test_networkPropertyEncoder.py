import os
import shutil
from unittest import TestCase
import random

from source.caching.CacheType import CacheType
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.data_model.receptor.receptor_sequence.SequenceMetadata import SequenceMetadata
from source.data_model.repertoire.Repertoire import Repertoire
from source.encodings.EncoderParams import EncoderParams
from source.encodings.network.NetworkPropertyEncoder import NetworkPropertyEncoder
from source.environment.Constants import Constants
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.environment.LabelConfiguration import LabelConfiguration
from source.util.PathBuilder import PathBuilder


class TestNetworkGlobalPropertyEncoder(TestCase):

    def setUp(self) -> None:
        os.environ[Constants.CACHE_TYPE] = CacheType.TEST.name

    def test_encode(self):
        path = EnvironmentSettings.root_path + "test/tmp/network/"

        PathBuilder.build(path)

        rep1 = Repertoire.build_from_sequence_objects(
            sequence_objects=[ReceptorSequence(x, metadata=SequenceMetadata(count=10)) for x in self.generate_random_sequences(100, 4, 20)],
            metadata={"l1": "test_1", "l2": 2}, path=path)

        rep2 = Repertoire.build_from_sequence_objects(
            sequence_objects=[ReceptorSequence(x, metadata=SequenceMetadata(count=10)) for x in self.generate_random_sequences(100, 4, 20)],
            metadata={"l1": "test_2", "l2": 3}, path=path)

        lc = LabelConfiguration()
        lc.add_label("l1", ["test_1", "test_2"])
        lc.add_label("l2", [0, 3])

        dataset = RepertoireDataset(repertoires=[rep1, rep2])

        encoder = NetworkPropertyEncoder.build_object(dataset, **{
            "max_edit_distance": 1
        })

        d1 = encoder.encode(dataset, EncoderParams(
            result_path=path + "1/",
            label_configuration=lc,
            batch_size=2
        ))

        shutil.rmtree(path)

    def generate_random_sequences(self, N, minimum_length=4, maximum_length=20, letters='ACDEFGHIKLMNPQRSTVWY'):

        string_list = []

        for i in range(N):
            random_integer = random.randrange(minimum_length, maximum_length)
            s = ''.join(random.choice(letters) for _ in range(random_integer))
            string_list.append(s)

        return string_list
