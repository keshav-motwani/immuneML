import hashlib
import math
import time
from multiprocessing.pool import Pool

import numpy as np

from source.analysis.network.NetworkConstructor import NetworkConstructor
from source.caching.CacheHandler import CacheHandler
from source.caching.CacheObjectType import CacheObjectType
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.encodings.EncoderParams import EncoderParams
from source.encodings.network.NetworkPropertyEncoder import NetworkPropertyEncoder
from source.logging.Logger import log


class NetworkPropertyRepertoireEncoder(NetworkPropertyEncoder):

    @log
    def _encode_new_dataset(self, dataset, params: EncoderParams):

        encoded_data = self._encode_data(dataset, params)

        encoded_dataset = RepertoireDataset(repertoires=dataset.repertoires,
                                            encoded_data=encoded_data,
                                            params=dataset.params,
                                            metadata_file=dataset.metadata_file)

        self.store(encoded_dataset, params)

        return encoded_dataset

    @log
    def _encode_examples(self, dataset, params: EncoderParams):

        arguments = [(repertoire, params) for repertoire in dataset.repertoires]

        with Pool(params["batch_size"]) as pool:
            chunksize = math.floor(dataset.get_example_count()/params["batch_size"]) + 1
            repertoires = pool.starmap(self.get_encoded_repertoire, arguments, chunksize=chunksize)

        encoded_repertoire_list, repertoire_names, labels = zip(*repertoires)

        encoded_labels = {k: [dic[k] for dic in labels] for k in labels[0]}

        return list(encoded_repertoire_list), list(repertoire_names), encoded_labels

    def get_encoded_repertoire(self, repertoire, params: EncoderParams):

        params["model"] = vars(self)

        return CacheHandler.memo_by_params((("encoding_model", params["model"]),
                                            ("labels", params["label_configuration"].get_labels_by_name()),
                                            ("repertoire_id", repertoire.identifier),
                                            ("repertoire_data",  hashlib.sha256(np.ascontiguousarray(repertoire.get_sequence_aas())).hexdigest())),
                                           lambda: self.encode_repertoire(repertoire, params), CacheObjectType.ENCODING_STEP)

    def encode_repertoire(self, repertoire, params: EncoderParams):

        a = time.time()

        frame_types = np.char.upper(repertoire.get_attribute("frame_types").astype(str))
        sequences = repertoire.get_attribute("sequence_aas").astype(str)
        sequences = sequences[frame_types == "IN"]

        graph = NetworkConstructor.construct_network(sequences, self.max_edit_distance)

        encoded = NetworkConstructor.compute_graph_properties(graph)["value"]

        label_config = params["label_configuration"]
        labels = dict()
        for label_name in label_config.get_labels_by_name():
            label = repertoire.metadata[label_name]
            labels[label_name] = label

        b = time.time()

        print(repertoire.data_filename, "took", b - a, "seconds")

        return encoded, repertoire.identifier, labels
