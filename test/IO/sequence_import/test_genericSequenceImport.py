import os
import shutil
from unittest import TestCase

from source.IO.sequence_import.GenericSequenceImport import GenericSequenceImport
from source.caching.CacheType import CacheType
from source.environment.Constants import Constants
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.util.PathBuilder import PathBuilder


class TestSequenceImport(TestCase):

    def setUp(self) -> None:
        os.environ[Constants.CACHE_TYPE] = CacheType.TEST.name

    def test_import_all_sequences(self):
        file_content = """Clone ID,Senior Author,TRAJ Gene,TRAV Gene,CDR3A AA Sequence,TRBV Gene,TRBD Gene,TRBJ Gene,CDR3B AA Sequence,Antigen Protein,Antigen Gene,Antigen Species,Antigen Peptide AA #,Epitope Peptide,MHC Class,HLA Restriction
1E6,Sewell,TRAJ12,TRAV12-3,CAMRGDSSYKLIF,TRBV12-4,TRBD2,TRBJ2-4,CASSLWEKLAKNIQYF,PPI,INS,Human,12-24,ALWGPDPAAA,MHC I,A*02:01
4.13,Nepom,TRAJ44,TRAV19,CALSENRGGTASKLTF,TRBV5-1,TRBD1,TRBJ1-1,CASSLVGGPSSEAFF,GAD,,Human,555-567,,MHC II,DRB1*04:01
5,Roep,TRAJ6,TRAV21,CAVKRTGGSYIPTF,TRBV11-2,TRBD1,TRBJ2-2,CASSSFWGSDTGELFF,Insulin B,,Human,9-23,,MHC II,DQ8
D222D 2,Mallone,TRAJ36*01,TRAV17*01,CAVTGANNLFF,TRBV19*01,TRBD1*01,TRBJ2-2*01,CASSIEGPTGELFF,Zinc Transporter 8,ZnT8,Human,185-194,AVAANIVLTV,MHC I,A*02:01
GSE.20D11,Nakayama,TRAJ4,TRAV12-3,CAILSGGYNKLIF,TRBV2,TRBD2,TRBJ2-5,CASSAETQYF,Insulin B,,Human,9-23,,MHC II,DQ8
GSE.6H9,Nakayama,TRAJ40,TRAV26-1,CIVRVDSGTYKYIF,TRBV7-2,TRBD2,TRBJ2-1,CASSLTAGLASTYNEQFF,Insulin B,,Human,9-23,,MHC II,DQ8/DQ8
iGRP 32,DiLorenzo,TRAJ48,TRAV12-1,CVVNILSNFGNEKLTF,TRBV20/OR9-2,TRBD1,TRBJ2-1,CSASRQGWVNEQFF,IGRP,,Human,265-273,,MHC I,A*02:01
MHB10.3,TBD,TRAJ27,TRAV4,CLVGDSLNTNAGKSTF,TRBV29-1,TRBD2,TRBJ2-2,CSVEDRNTGELFF,Insulin B,,Human,11-30,,,DRB1*03:01
PM1#11,TBD,TRAJ54,TRAV35,CAGHSIIQGAQKLVF,TRBV5-1,TRBD2,TRBJ2-1,CASGRSSYNEQFF,GAD,,Human,339-352,,,DRB1*03:01
R164,Nepom,TRAJ56,TRAV19,CALSEEGGGANSKLTF,TRBV5-1,TRBD2,TRBJ1-6,CASSLAGGANSPLHF,GAD,,Human,555-567,,,DRB1*04:01
SD32.5,Boehm,TRAJ23,TRAV26-1,CIVRVSSAYYNQGGKLIF,TRBV27,TRBD2,TRBJ2-3,CASSPRANTDTQYF,Insulin A,,Human,5-21,,,DRB1*04:01
SD52.c1,Boehm,TRAJ27,TRAV4,CLVGDSLNTNAGKSTF,TRBV27,TRBD1,TRBJ1-5,CASSWSSIGNQPQHF,PPI,INS,Human,C18-A1,,,DRB1*04:01
T1D#10 C8,TBD,TRAJ26,TRAV12-3,CATAYGQNFVF,TRBV4-1,TRBD2,TRBJ2-2,CASSRGGGNTGELFF,Insulin B,,Human,9-23,,MHC II,DQ8
T1D#3 C8,TBD,TRAJ23,TRAV17,CATDAGYNQGGKLIF,TRBV5-1,TRBD2,TRBJ1-3,CASSAGNTIYF,Insulin B,,Human,9-23,,MHC II,DQ8"""

        path = EnvironmentSettings.root_path + "test/tmp/sequence_import/"
        PathBuilder.build(path)

        with open(path + "seqs.tsv", "w") as file:
            file.writelines(file_content)

        sequences = GenericSequenceImport.import_items(path + "seqs.tsv",
                                                       column_mapping={"CDR3B AA Sequence": "amino_acid",
                                                                       "TRBV Gene": "v_gene",
                                                                       "TRBJ Gene": "j_gene"},
                                                       separator=",")

        self.assertEqual(14, len(sequences))
        self.assertEqual("ASSLWEKLAKNIQY", sequences[0].get_sequence())
        self.assertEqual("TRBJ2-2*01", sequences[3].metadata.j_gene)
        self.assertEqual("TRBV19*01", sequences[3].metadata.v_gene)

        shutil.rmtree(path)
