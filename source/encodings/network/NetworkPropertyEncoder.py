import abc

import numpy as np
import pandas as pd
from scipy import sparse

from source.IO.dataset_export.PickleExporter import PickleExporter
from source.analysis.network.NetworkConstructor import NetworkConstructor
from source.caching.CacheHandler import CacheHandler
from source.data_model.encoded_data.EncodedData import EncodedData
from source.encodings.DatasetEncoder import DatasetEncoder
from source.encodings.EncoderParams import EncoderParams
from source.util.ReflectionHandler import ReflectionHandler


class NetworkPropertyEncoder(DatasetEncoder):
    """
    The EvennessProfileEncoder class encodes a repertoire based on the clonal frequency distribution. The evenness for
    a given repertoire is defined as follows:

    .. math::

        ^{\alpha} \mathrm{E}(\mathrm{f})=\frac{\left(\sum_{\mathrm{i}=1}^{\mathrm{n}} \mathrm{f}_{\mathrm{i}}^{\alpha}\right)^{\frac{1}{1-\alpha}}}{\mathrm{n}}

    That is, it is the exponential of Renyi entropy at a given alpha divided by the species richness, or number of unique
    sequences.

    See Greiff and colleagues' publication "A bioinformatic framework for immune repertoire diversity profiling enables
    detection of immunological status" in Genome Medicine 2015 for more details.

    Attributes:

        min_alpha (float): minimum alpha value to use

        max_alpha (float): maximum alpha value to use

        dimension (int): dimension of output evenness profile vector, or the number of alpha values to linearly space
            between min_alpha and max_alpha

    Specification:

    .. indent with spaces
    .. code-block:: yaml

            my_evenness_profile:
                EvennessProfile:
                    min_alpha: 0
                    max_alpha: 10
                    dimension: 51


    """

    STEP_ENCODED = "encoded"
    STEP_VECTORIZED = "vectorized"

    dataset_mapping = {
        "RepertoireDataset": "NetworkPropertyRepertoireEncoder",
    }

    def __init__(self, max_edit_distance: int = 1, name: str = None):
        self.max_edit_distance = max_edit_distance
        self.name = name

    @staticmethod
    def _prepare_parameters(max_edit_distance: int, name: str = None):

        return {
            "max_edit_distance": max_edit_distance,
            "name": name
        }

    @staticmethod
    def build_object(dataset=None, **params):
        try:
            prepared_params = NetworkPropertyEncoder._prepare_parameters(**params)
            encoder = ReflectionHandler.get_class_by_name(NetworkPropertyEncoder.dataset_mapping[dataset.__class__.__name__],
                                                          "network/")(**prepared_params)
        except ValueError:
            raise ValueError("{} is not defined for dataset of type {}.".format(NetworkPropertyEncoder.__name__, dataset.__class__.__name__))
        return encoder

    def encode(self, dataset, params: EncoderParams):

        encoded_dataset = CacheHandler.memo_by_params(self._prepare_caching_params(dataset, params),
                                                      lambda: self._encode_new_dataset(dataset, params))

        return encoded_dataset

    def _prepare_caching_params(self, dataset, params: EncoderParams, step: str = ""):
        return (("example_identifiers", tuple(dataset.get_example_ids())),
                ("dataset_type", dataset.__class__.__name__),
                ("labels", tuple(params["label_configuration"].get_labels_by_name())),
                ("encoding", NetworkPropertyEncoder.__name__),
                ("learn_model", params["learn_model"]),
                ("step", step),
                ("encoding_params", tuple(vars(self).items())))

    def _encode_data(self, dataset, params: EncoderParams) -> EncodedData:

        encoded_example_list, example_ids, encoded_labels = CacheHandler.memo_by_params(
            self._prepare_caching_params(dataset, params, NetworkPropertyEncoder.STEP_ENCODED),
            lambda: self._encode_examples(dataset, params))

        vectorized_examples = CacheHandler.memo_by_params(
            self._prepare_caching_params(dataset, params, NetworkPropertyEncoder.STEP_VECTORIZED),
            lambda: self._vectorize_encoded(examples=encoded_example_list))

        feature_names = NetworkConstructor.GRAPH_PROPERTIES

        feature_annotations = pd.DataFrame({"feature": feature_names})

        encoded_data = EncodedData(examples=vectorized_examples,
                                   labels=encoded_labels,
                                   feature_names=feature_names,
                                   example_ids=example_ids,
                                   feature_annotations=feature_annotations,
                                   encoding=NetworkPropertyEncoder.__name__)

        return encoded_data

    @abc.abstractmethod
    def _encode_new_dataset(self, dataset, params: EncoderParams):
        pass

    @abc.abstractmethod
    def _encode_examples(self, dataset, params: EncoderParams):
        pass

    def store(self, encoded_dataset, params: EncoderParams):
        PickleExporter.export(encoded_dataset, params["result_path"])

    def _vectorize_encoded(self, examples: list):

        vectorized_examples = sparse.csr_matrix(np.vstack(examples))

        return vectorized_examples
