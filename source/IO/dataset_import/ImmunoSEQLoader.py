import os

from source.IO.dataset_import.DataLoader import DataLoader
from source.IO.dataset_import.PickleLoader import PickleLoader
from source.data_model.dataset.Dataset import Dataset
from source.util.FilenameHandler import FilenameHandler


class ImmunoSEQLoader(DataLoader):

    COLUMNS = ['nucleotide', 'aminoAcid', 'count (tempates/reads)', 'vMaxResolved', 'jMaxResolved', 'vFamilyName', 'jFamilyName', 'sequenceStatus']

    @staticmethod
    def load(path, params: dict = None) -> Dataset:

        file_path = path + FilenameHandler.get_dataset_name(ImmunoSEQLoader.__name__)

        if os.path.isfile(file_path):
            dataset = PickleLoader.load(file_path)
        else:
            dataset = None
        #     dataset = ImmunoSEQLoader._process_dataset(path, params)

        return dataset

    # @staticmethod
    # def _process_dataset(path, params: dict = None):
    #
    #     filenames = ImmunoSEQLoader._get_filenames(path, params)
    #     PathBuilder.build(params["result_path"])
    #     repertoire_filenames, sample_params = ImmunoSEQLoader._process_repertoires(filenames, params)
    #
    #     dataset = Dataset(filenames=repertoire_filenames, params=sample_params)
    #
    #     PickleExporter.export(dataset, path, FilenameHandler.get_dataset_name(ImmunoSEQLoader.__name__))
    #
    #     return dataset
    #
    # @staticmethod
    # def _process_repertoire(filename, params):
    #     basename = os.path.splitext(os.path.basename(filename))[0]
    #     filepath = params["result_path"] + basename + ".pickle"
    #
    #     if not os.path.isfile(filepath):
    #         ImmunoSEQLoader._load_repertoire(filename, params, filepath)
    #
    #     return filepath

    # @staticmethod
    # def _load_repertoire(filename, params, result_file_path):
    #     repertoire = Repertoire()
    #
    #     # add metadata
    #     sample = SampleParser.parse(filename, params["sample_name_parser"], params["sample_parser_params"])
    #     repertoire.metadata = RepertoireMetadata(sample=sample)
    #
    #     # add sequences
    #     sequences = ImmunoSEQLoader._preprocess_sequences()
    #     repertoire.sequences = sequences
    #
    #     # store as pickle
    #     with open(result_file_path, "wb") as file:
    #         pickle.dump(repertoire, file)
    #
    #
    # @staticmethod
    # def _process_repertoires(filenames: list, params: dict):
    #     process_count = ParallelismManager.assign_cores_to_job("load_experimental_data")
    #
    #     with Pool(process_count) as pool:
    #         results = pool.starmap(ImmunoSEQLoader._process_repertoire, product(filenames, params))
    #
    #     outputs = [result[0] for result in results]
    #
    #     return outputs
    #
    # @staticmethod
    # def _get_filenames(path, params: dict):
    #     ext = params["file_type"]  # should be tsv, csv...
    #     if "file_names" in params and params["file_names"] is not None:
    #         filenames = params["file_names"]
    #     else:
    #         filenames = sorted(glob.glob(path + "*." + ext))
    #
    #     if "file_indices" is not None:
    #         filenames = [filenames[i] for i in params["file_indices"]]
    #
    #     return filenames