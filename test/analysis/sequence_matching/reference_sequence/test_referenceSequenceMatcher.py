import shutil
from unittest import TestCase

from source.analysis.sequence_matching.reference_sequence.ReferenceSequenceMatcher import ReferenceSequenceMatcher
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.data_model.receptor.receptor_sequence.SequenceMetadata import SequenceMetadata

from source.data_model.repertoire.Repertoire import Repertoire
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.util.PathBuilder import PathBuilder


class TestReferenceSequenceMatcher(TestCase):

    def test_match(self):

        path = EnvironmentSettings.root_path + "test/tmp/seqmatch/"
        PathBuilder.build(path)

        repertoire = Repertoire.build_from_sequence_objects(sequence_objects=[
            ReceptorSequence(amino_acid_sequence="AAAAAA", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J2"), identifier="3"),
            ReceptorSequence(amino_acid_sequence="CCCCCC", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J2"), identifier="4"),
            ReceptorSequence(amino_acid_sequence="AAAACC", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J2"), identifier="5"),
            ReceptorSequence(amino_acid_sequence="TADQVF", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J3"), identifier="6")],
            metadata={"CD": True}, path=path)

        dataset = RepertoireDataset(repertoires=[repertoire])
        sequences = [ReceptorSequence("AAAACA", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J2"), identifier="1"),
                     ReceptorSequence("TADQV", metadata=SequenceMetadata(chain="A", v_gene="V1", j_gene="J3"), identifier="2")]

        matcher = ReferenceSequenceMatcher
        result = matcher.match(dataset, sequences, False, ["v_gene", "j_gene"], 1, batch_size=8)

        self.assertTrue(len(result.repertoires) == 1)
        self.assertTrue(len(result.repertoires[0].reference_sequences) == 2) # number of reference sequences
        self.assertTrue(len(result.repertoires[0].reference_sequences[0].matching_query_sequences) == 2) # number of matching sequences for 1st reference sequence
        self.assertTrue(len(result.repertoires[0].reference_sequences[1].matching_query_sequences) == 1) # number of matching sequences for 2nd reference sequence

        shutil.rmtree(path)
