from typing import List

from source.IO.dataset_export.DataExporter import DataExporter
from source.dsl.symbol_table.SymbolTable import SymbolTable
from source.dsl.symbol_table.SymbolType import SymbolType
from source.util.ParameterValidator import ParameterValidator
from source.util.ReflectionHandler import ReflectionHandler
from source.workflows.instructions.dataset_generation.DatasetGenerationInstruction import DatasetGenerationInstruction


class DatasetGenerationParser:
    """
    Specification of instruction with a random datasets:

    definitions:
      datasets:
        my_generated_dataset: # a dataset to be exported in the given format
          format: RandomRepertoireDataset
          params:
            result_path: generated_dataset/
            repertoire_count: 100
            sequence_count_probabilities:
              100: 0.5
              120: 0.5
            sequence_length_probabilities:
              12: 0.333
              13: 0.333
              14: 0.333
            labels:
              immune_event_1:
                yes: 0.5
                no: 0.5
    instructions:
      my_instruction1: # instruction name
        type: DatasetGeneration
        datasets: # list of datasets to export
          - my_generated_dataset
        formats: # list of formats to export the datasets to
          - AIRR
          - Pickle
    """

    VALID_KEYS = ["type", "datasets", "formats"]

    def parse(self, key: str, instruction: dict, symbol_table: SymbolTable) -> DatasetGenerationInstruction:
        location = "DatasetGenerationParser"
        ParameterValidator.assert_keys(list(instruction.keys()), DatasetGenerationParser.VALID_KEYS, location, key)
        valid_formats = DatasetGenerationParser.get_valid_formats()
        ParameterValidator.assert_all_in_valid_list(instruction["formats"], valid_formats, location, "formats")
        ParameterValidator.assert_all_in_valid_list(instruction["datasets"], symbol_table.get_keys_by_type(SymbolType.DATASET), location,
                                                    "datasets")

        return DatasetGenerationInstruction(datasets=[symbol_table.get(dataset_key) for dataset_key in instruction["datasets"]],
                                            exporters=[ReflectionHandler.get_class_by_name(f"{key}Exporter", "dataset_export/")
                                                       for key in instruction["formats"]],
                                            name=key)

    @staticmethod
    def get_valid_formats() -> List[str]:
        class_path = "dataset_export/"
        classes = ReflectionHandler.get_classes_by_partial_name("Exporter", class_path)
        valid_values = [cls.__name__[:-8] for cls in ReflectionHandler.all_nonabstract_subclasses(DataExporter)]
        return valid_values
