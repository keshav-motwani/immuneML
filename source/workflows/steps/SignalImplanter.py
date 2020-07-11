import copy
import os
from typing import List

import pandas as pd

from source.IO.dataset_import.PickleImport import PickleImport
from source.data_model.dataset.Dataset import Dataset
from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.data_model.repertoire.Repertoire import Repertoire
from source.simulation.SimulationState import SimulationState
from source.util.FilenameHandler import FilenameHandler
from source.util.PathBuilder import PathBuilder
from source.workflows.steps.Step import Step


class SignalImplanter(Step):

    DATASET_NAME = "simulated_dataset"

    @staticmethod
    def run(input_params: SimulationState = None):
        path = input_params.result_path + FilenameHandler.get_dataset_name(SignalImplanter.__name__)

        if os.path.isfile(path):
            dataset = PickleImport.import_dataset({"path": path}, SignalImplanter.DATASET_NAME)
        else:
            dataset = SignalImplanter._implant_signals(input_params)

        return dataset

    @staticmethod
    def _implant_signals(input_params: SimulationState = None) -> Dataset:

        PathBuilder.build(input_params.result_path)
        PathBuilder.build(input_params.result_path + "repertoires/")

        processed_repertoires = []
        simulation_limits = SignalImplanter._prepare_simulation_limits(input_params.simulation.implantings,
                                                                       input_params.dataset.get_example_count())
        simulation_index = 0

        implanting_metadata = {f"signal_{signal.id}": [] for signal in input_params.signals}

        for index, repertoire in enumerate(input_params.dataset.get_data(input_params.batch_size)):

            if simulation_index <= len(simulation_limits) - 1 and index >= simulation_limits[simulation_index]:
                simulation_index += 1

            processed_repertoire = SignalImplanter._process_repertoire(index, repertoire, simulation_index, simulation_limits, input_params)
            processed_repertoires.append(processed_repertoire)

        processed_dataset = RepertoireDataset(repertoires=processed_repertoires, params=input_params.dataset.params, name=input_params.dataset.name,
                                              metadata_file=SignalImplanter._create_metadata_file(processed_repertoires, input_params))
        return processed_dataset

    @staticmethod
    def _create_metadata_file(processed_repertoires: List[Repertoire], input_params) -> str:

        path = input_params.result_path + "metadata.csv"

        new_df = pd.DataFrame([repertoire.metadata for repertoire in processed_repertoires])
        new_df.drop('field_list', axis=1, inplace=True)
        new_df["filename"] = [repertoire.data_filename for repertoire in processed_repertoires]
        new_df.to_csv(path, index=False)

        return path

    @staticmethod
    def _process_repertoire(index, repertoire, simulation_index, simulation_limits, input_params) -> Repertoire:
        if simulation_index < len(simulation_limits):
            return SignalImplanter._implant_in_repertoire(index, repertoire, simulation_index, input_params)
        else:
            return SignalImplanter._copy_repertoire(index, repertoire, input_params)

    @staticmethod
    def _copy_repertoire(index: int, repertoire: Repertoire, input_params: SimulationState) -> Repertoire:
        new_repertoire = Repertoire.build_from_sequence_objects(repertoire.sequences, input_params.result_path + "repertoires/", repertoire.metadata)

        for signal in input_params.signals:
            new_repertoire.metadata[f"signal_{signal.id}"] = False

        return new_repertoire

    @staticmethod
    def _implant_in_repertoire(index, repertoire, simulation_index, input_params) -> Repertoire:
        new_repertoire = copy.deepcopy(repertoire)
        for signal in input_params.simulation.implantings[simulation_index].signals:
            new_repertoire = signal.implant_to_repertoire(repertoire=new_repertoire,
                                                          repertoire_implanting_rate=
                                                          input_params.simulation.implantings[simulation_index].repertoire_implanting_rate,
                                                          path=input_params.result_path + "repertoires/")

        for signal in input_params.simulation.implantings[simulation_index].signals:
            new_repertoire.metadata[f"signal_{signal.id}"] = True
        for signal in input_params.signals:
            if signal not in input_params.simulation.implantings[simulation_index].signals:
                new_repertoire.metadata[f"signal_{signal.id}"] = False

        return new_repertoire

    @staticmethod
    def _prepare_simulation_limits(simulation: list, repertoire_count: int) -> list:
        limits = [int(item.dataset_implanting_rate * repertoire_count) for item in simulation]
        limits = [sum(limits[:i+1]) for i in range(len(limits))]
        return limits
