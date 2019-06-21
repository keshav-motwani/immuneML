import copy
import os

import numpy as np
import pandas as pd

from source.analysis.criteria_matches.CriteriaMatcher import CriteriaMatcher
from source.data_model.dataset.Dataset import Dataset
from source.preprocessing.Preprocessor import Preprocessor


class MetadataFilter(Preprocessor):
    """
    For use in filtering out repertoires from the dataset based on information stored in the metadata_file.

    Requires params with only one key: "criteria" which is to be specified based on the format in CriteriaMatcher

    For example:

    params = {
        "criteria": {
            "type": OperationType.GREATER_THAN,
            "value": {
                "type": DataType.COLUMN,
                "name": "key2"
            },
            "threshold": 1
        }
    }

    This filter includes only repertoires with values greater than 1 in the "key2" column of the metadata_file
    """

    @staticmethod
    def process(dataset: Dataset, params: dict) -> Dataset:
        processed_dataset = copy.deepcopy(dataset)
        original_filenames = processed_dataset.get_filenames()
        indices = MetadataFilter.get_matching_indices(processed_dataset, params["criteria"])
        processed_dataset.set_filenames([original_filenames[i] for i in indices])
        processed_dataset.metadata_file = MetadataFilter.build_new_metadata(dataset, indices)
        return processed_dataset

    @staticmethod
    def get_matching_indices(dataset: Dataset, criteria):
        metadata = pd.DataFrame(dataset.get_metadata(None))
        matches = CriteriaMatcher().match(criteria, metadata)
        indices = np.where(matches)[0]
        return indices

    @staticmethod
    def build_new_metadata(dataset, indices_to_keep):
        if dataset.metadata_file:
            df = pd.read_csv(dataset.metadata_file, index_col=0).iloc[indices_to_keep, :]
            path = os.path.dirname(os.path.abspath(dataset.metadata_file)) + "/{}_metadata_filtered.csv"\
                .format(os.path.splitext(os.path.basename(dataset.metadata_file))[0])
            df.to_csv(path)
        else:
            path = None
        return path