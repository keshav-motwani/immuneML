import shutil
from unittest import TestCase

import numpy as np
import pandas as pd

from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.encodings.EncoderParams import EncoderParams
from source.encodings.reference_encoding.MatchedReferenceSequenceEncoder import MatchedReferenceSequenceEncoder
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.environment.LabelConfiguration import LabelConfiguration
from source.util.RepertoireBuilder import RepertoireBuilder


class TestMatchedReferenceSequenceEncoder(TestCase):
    def test__encode_new_dataset(self):
        path = EnvironmentSettings.root_path + "test/tmp/matched_ref_encoder/"
        repertoires, metadata = RepertoireBuilder.build([["AAAA", "AACA"], ["TTTA", "AAAA"]], path, {"default": np.array([1, 2])})
        dataset = RepertoireDataset(repertoires=repertoires)

        label_config = LabelConfiguration()
        label_config.add_label("default", [1, 2])

        file_content = """complex.id	Gene	CDR3	V	J	Species	MHC A	MHC B	MHC class	Epitope	Epitope gene	Epitope species	Reference	Method	Meta	CDR3fix	Score
        100a	TRA	AAAA	TRAV12	TRAJ1	HomoSapiens	HLA-A*11:01	B2M	MHCI	AVFDRKSDAK	EBNA4	EBV	https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#	{"frequency": "1/11684", "identification": "dextramer-sort", "sequencing": "rna-seq", "singlecell": "yes", "verification": ""}	{"cell.subset": "", "clone.id": "", "donor.MHC": "", "donor.MHC.method": "", "epitope.id": "", "replica.id": "", "samples.found": 1, "structure.id": "", "studies.found": 1, "study.id": "", "subject.cohort": "", "subject.id": "1", "tissue": ""}	{"cdr3": "CASSPPRVYSNGAGLAGVGWRNEQFF", "cdr3_old": "CASSPPRVYSNGAGLAGVGWRNEQFF", "fixNeeded": false, "good": true, "jCanonical": true, "jFixType": "NoFixNeeded", "jId": "TRBJ2-1*01", "jStart": 21, "vCanonical": true, "vEnd": 4, "vFixType": "NoFixNeeded", "vId": "TRBV5-4*01"}	0
        100a	TRA	BBBB	TRAV12	TRAJ1	HomoSapiens	HLA-A*11:01	B2M	MHCI	AVFDRKSDAK	EBNA4	EBV	https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#	{"frequency": "1/11684", "identification": "dextramer-sort", "sequencing": "rna-seq", "singlecell": "yes", "verification": ""}	{"cell.subset": "", "clone.id": "", "donor.MHC": "", "donor.MHC.method": "", "epitope.id": "", "replica.id": "", "samples.found": 1, "structure.id": "", "studies.found": 1, "study.id": "", "subject.cohort": "", "subject.id": "1", "tissue": ""}	{"cdr3": "CASSPPRVYSNGAGLAGVGWRNEQFF", "cdr3_old": "CASSPPRVYSNGAGLAGVGWRNEQFF", "fixNeeded": false, "good": true, "jCanonical": true, "jFixType": "NoFixNeeded", "jId": "TRBJ2-1*01", "jStart": 21, "vCanonical": true, "vEnd": 4, "vFixType": "NoFixNeeded", "vId": "TRBV5-4*01"}	0
        100a	TRA	CCCC	TRAV12	TRAJ1	HomoSapiens	HLA-A*11:01	B2M	MHCI	AVFDRKSDAK	EBNA4	EBV	https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#	{"frequency": "1/11684", "identification": "dextramer-sort", "sequencing": "rna-seq", "singlecell": "yes", "verification": ""}	{"cell.subset": "", "clone.id": "", "donor.MHC": "", "donor.MHC.method": "", "epitope.id": "", "replica.id": "", "samples.found": 1, "structure.id": "", "studies.found": 1, "study.id": "", "subject.cohort": "", "subject.id": "1", "tissue": ""}	{"cdr3": "CASSPPRVYSNGAGLAGVGWRNEQFF", "cdr3_old": "CASSPPRVYSNGAGLAGVGWRNEQFF", "fixNeeded": false, "good": true, "jCanonical": true, "jFixType": "NoFixNeeded", "jId": "TRBJ2-1*01", "jStart": 21, "vCanonical": true, "vEnd": 4, "vFixType": "NoFixNeeded", "vId": "TRBV5-4*01"}	0
        """

        with open(path + "refs.tsv", "w") as file:
            file.writelines(file_content)

        references = {"path": path + "refs.tsv", "format": "VDJdb"}

        encoder = MatchedReferenceSequenceEncoder.build_object(dataset, **{
                "reference_sequences": references,
                "max_edit_distance": 1,
                "same_length_sequence": True,
                "metadata_attrs_to_match": [],
                "summary": "PCT_UNIQUE_READS_WITH_MATCH"
            })

        encoded = encoder.encode(dataset, EncoderParams(
            result_path=path,
            label_configuration=label_config,
            model={},
            filename="dataset.csv"
        ))

        self.assertEqual((2, 3), encoded.encoded_data.examples.shape)
        self.assertTrue(isinstance(encoder, MatchedReferenceSequenceEncoder))
        self.assertTrue(encoded.encoded_data.examples[0, 0] == 1)
        self.assertTrue(encoded.encoded_data.examples[1, 0] == 0.5)
        self.assertTrue(encoded.encoded_data.labels["default"] == [1, 2])
        self.assertTrue(encoded.encoded_data.feature_names == ["AAAA()", "BBBB()", "CCCC()"])
        self.assertTrue(isinstance(encoded.encoded_data.feature_annotations, pd.DataFrame))

        shutil.rmtree(path)
