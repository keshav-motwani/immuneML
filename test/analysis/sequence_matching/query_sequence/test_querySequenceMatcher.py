import os
import shutil
from unittest import TestCase

from source.analysis.sequence_matching.query_sequence.QuerySequenceMatcher import QuerySequenceMatcher
from source.caching.CacheType import CacheType
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.data_model.receptor.receptor_sequence.SequenceMetadata import SequenceMetadata

from source.data_model.repertoire.Repertoire import Repertoire
from source.environment.Constants import Constants
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.util.PathBuilder import PathBuilder


class TestQuerySequenceMatcher(TestCase):

    def setUp(self) -> None:
        os.environ[Constants.CACHE_TYPE] = CacheType.TEST.name

    def test_match(self):

        path = EnvironmentSettings.root_path + "test/tmp/seqmatch/"
        PathBuilder.build(path)

        repertoire = Repertoire.build_from_sequence_objects(sequence_objects=
                                                                    [ReceptorSequence(amino_acid_sequence="AAAAAA", identifier="1",
                                                                                      metadata=SequenceMetadata(chain="A", count=3)),
                                                                     ReceptorSequence(amino_acid_sequence="CCCCCC", identifier="2",
                                                                                      metadata=SequenceMetadata(chain="A", count=2)),
                                                                     ReceptorSequence(amino_acid_sequence="AAAACC", identifier="3",
                                                                                      metadata=SequenceMetadata(chain="A", count=1)),
                                                                     ReceptorSequence(amino_acid_sequence="TADQVF", identifier="4",
                                                                                      metadata=SequenceMetadata(chain="A", count=4))],
                                                            metadata={"CD": True}, path=path)

        sequences = [ReceptorSequence("AAAACA", metadata=SequenceMetadata(chain="A")),
                     ReceptorSequence("TADQV", metadata=SequenceMetadata(chain="A"))]

        dataset = RepertoireDataset(repertoires=[repertoire])

        matcher = QuerySequenceMatcher()
        result = matcher.match(dataset, sequences, False, ["v_gene", "j_gene", "chain"], 1, batch_size=1)

        self.assertEqual(1, len(result.repertoires[0].query_sequences[3].matching_reference_sequences)) # number of matching sequences for the 4th sequence in the 1st repertoire
        self.assertTrue(result.repertoires[0].metadata["CD"])
        self.assertEqual(1, len(result.repertoires))

        result = result.repertoires[0]

        self.assertTrue(result.query_sequences is not None)
        self.assertTrue(result.identifier is not None)
        self.assertTrue(result.index == 0)

        self.assertEqual(4, len(result.query_sequences)) # number of sequences in the repertoire
        self.assertEqual(1, len(result.query_sequences[0].matching_reference_sequences)) # number of matching sequences for the 1st sequence in the repertoire
        self.assertEqual(0, len(result.query_sequences[1].matching_reference_sequences))
        self.assertEqual(1, len(result.query_sequences[2].matching_reference_sequences))
        self.assertEqual(1, len(result.query_sequences[3].matching_reference_sequences))

        self.assertEqual(3, len([r for r in result.query_sequences if len(r.matching_reference_sequences) > 0])) # number of sequences in repertoire with at least one match
        self.assertTrue(result.metadata["CD"])
        self.assertEqual(0.8, result.pct_total_reads_with_match)

        shutil.rmtree(path)
