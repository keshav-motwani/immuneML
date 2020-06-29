import pandas as pd

from source.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from source.data_model.receptor.receptor_sequence.SequenceFrameType import SequenceFrameType
from source.data_model.receptor.receptor_sequence.SequenceMetadata import SequenceMetadata
from source.environment.Constants import Constants


class GenericSequenceImport:

    @staticmethod
    def import_items(path: str, params: dict):
        return GenericSequenceImport._read_sequences(path, params)

    @staticmethod
    def _read_sequences(filepath, params):

        usecols = None if params["additional_columns"] == "*" else list(params["column_mapping"].values()) + params["additional_columns"]
        separator = params.get("separator", "\t")

        try:
            df = pd.read_csv(filepath, sep=separator, iterator=False, usecols=usecols)
        except:
            df = pd.read_csv(filepath, sep=separator, iterator=False, usecols=usecols, encoding="latin1")

        df = df.rename(columns={j: i for i, j in params["column_mapping"].items()})

        if params.get("region_definition") is "IMGT":
            if "amino_acid" in df.columns:
                df['amino_acid'] = df["amino_acid"].str[1:-1]
            if "nucleotide" in df.columns:
                df['nucleotide'] = df["nucleotide"].str[3:-3]

        df = df.replace(["unresolved", "no data", "na", "unknown", "null", "nan", np.nan], Constants.UNKNOWN)

        return df.apply(GenericSequenceImport.create_sequence_from_row, axis=1, args=(params,)).values

    @staticmethod
    def create_sequence_from_row(row, params) -> ReceptorSequence:

        metadata = SequenceMetadata(v_subgroup=row.get("v_subgroup", Constants.UNKNOWN),
                                    v_gene=row.get("v_gene", Constants.UNKNOWN),
                                    v_allele=row.get("v_allele", Constants.UNKNOWN),
                                    j_subgroup=row.get("j_subgroup", Constants.UNKNOWN),
                                    j_gene=row.get("j_gene", Constants.UNKNOWN),
                                    j_allele=row.get("j_allele", Constants.UNKNOWN),
                                    chain=row.get("chain", "TRB"),
                                    count=int(row.get("count", "0")) if str(row.get("count", "0")).isdigit() else 0,
                                    frame_type=row.get("frame_type", SequenceFrameType.IN.value),
                                    region_type=row.get("region_type", "CDR3"))

        sequence = ReceptorSequence(amino_acid_sequence=row.get("amino_acid", None),
                                    nucleotide_sequence=row.get("nucleotide", None),
                                    metadata=metadata)

        for column in row.keys():
            if params["additional_columns"] == "*" or column in params["additional_columns"]:
                metadata.custom_params[column] = row[column]

        return sequence
