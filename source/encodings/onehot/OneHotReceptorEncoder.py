import numpy as np

from source.data_model.dataset.ReceptorDataset import ReceptorDataset
from source.data_model.encoded_data.EncodedData import EncodedData
from source.encodings.EncoderParams import EncoderParams
from source.encodings.onehot.OneHotEncoder import OneHotEncoder


class OneHotReceptorEncoder(OneHotEncoder):
    """
    One-hot encoded repertoire data is represented in a matrix with dimensions:
        [receptors, chains, one_hot_characters, sequence_lengths]

    when use_positional_info is true, the last 3 indices in one_hot_characters represents the positional information:
        - start position (high when close to start)
        - middle position (high in the middle of the sequence)
        - end position (high when close to end)
    """

    def _encode_new_dataset(self, dataset, params: EncoderParams):
        encoded_data = self._encode_data(dataset, params)
        if self.use_positional_info:
            encoded_data.examples = np.swapaxes(encoded_data.examples, 2, 3)

        encoded_dataset = ReceptorDataset(filenames=dataset.get_filenames(),
                                          encoded_data=encoded_data,
                                          params=dataset.params)

        self.store(encoded_dataset, params)

        return encoded_dataset

    def _encode_data(self, dataset: ReceptorDataset, params: EncoderParams):
        receptor_objs = [receptor for receptor in dataset.get_data(params["batch_size"])]
        sequences = [[getattr(obj, chain).get_sequence() for chain in obj.get_chains()] for obj in receptor_objs]
        first_chain_seqs, second_chain_seqs = zip(*sequences)

        max_seq_len = max(max([len(seq) for seq in first_chain_seqs]), max([len(seq) for seq in second_chain_seqs]))

        example_ids = dataset.get_example_ids()
        labels = self._get_labels(receptor_objs, params)

        examples_first_chain = self._encode_sequence_list(first_chain_seqs, pad_n_sequences=len(receptor_objs),
                                                          pad_sequence_len=max_seq_len)
        examples_second_chain = self._encode_sequence_list(second_chain_seqs, pad_n_sequences=len(receptor_objs),
                                                           pad_sequence_len=max_seq_len)

        examples = np.stack((examples_first_chain, examples_second_chain), axis=1)

        encoded_data = EncodedData(examples=examples,
                                   labels=labels,
                                   example_ids=example_ids,
                                   encoding=OneHotEncoder.__name__)

        return encoded_data

    def _get_labels(self, receptor_objs, params: EncoderParams):
        label_names = params["label_configuration"].get_labels_by_name()
        labels = {name: [None] * len(receptor_objs) for name in label_names}

        for idx, receptor in enumerate(receptor_objs):
            for name in label_names:
                labels[name][idx] = receptor.metadata[name]

        return labels
