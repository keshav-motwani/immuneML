import argparse
import os
import shutil
import sys; sys.path.extend([os.path.abspath(__file__ + "/../../../")])

from source.caching.CacheType import CacheType
from source.dsl.ImmuneMLParser import ImmuneMLParser
from source.dsl.semantic_model.SemanticModel import SemanticModel
from source.dsl.symbol_table.SymbolType import SymbolType
from source.environment.Constants import Constants
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.util.PathBuilder import PathBuilder
from source.util.ReflectionHandler import ReflectionHandler


class ImmuneMLApp:

    def __init__(self, specification_path: str, result_path: str):
        self._specification_path = specification_path
        self._result_path = os.path.relpath(result_path) + "/"
        self._cache_path = f"{self._result_path}cache/"

    def set_cache(self):
        os.environ[Constants.CACHE_TYPE] = CacheType.PRODUCTION.value
        EnvironmentSettings.set_cache_path(self._cache_path)

    def clear_cache(self):
        shutil.rmtree(self._cache_path, ignore_errors=True)
        EnvironmentSettings.reset_cache_path()
        del os.environ[Constants.CACHE_TYPE]

    def set_logging(self):
        if "ImmuneML_with_Galaxy" in os.environ and os.environ["ImmuneML_with_Galaxy"]:
            sys.stderr = open(self._result_path + "log.txt", 'w')

    def run(self):

        if self._result_path is not None:
            PathBuilder.build(self._result_path, warn_if_exists=True)

        self.set_logging()
        self.set_cache()

        symbol_table, self._specification_path = ImmuneMLParser.parse_yaml_file(self._specification_path,
                                                                                self._result_path)

        print("ImmuneML: starting the analysis...")

        instructions = symbol_table.get_by_type(SymbolType.INSTRUCTION)
        output = symbol_table.get("output")
        model = SemanticModel([instruction.item for instruction in instructions], self._result_path, output)
        result = model.run()

        self.clear_cache()
        return result


def run_immuneML(namespace: argparse.Namespace):
    if namespace.tool is None:
        app = ImmuneMLApp(namespace.yaml_path, namespace.output_dir)
    else:
        app_cls = ReflectionHandler.get_class_by_name(namespace.tool, "api/")
        app = app_cls(**vars(namespace))
    app.run()


def main():
    parser = argparse.ArgumentParser(description="immuneML command line tool")
    parser.add_argument("yaml_path", help="Path to specification YAML file. Always used to define the analysis.")
    parser.add_argument("output_dir", help="Output directory path.")
    parser.add_argument("--inputs", type=str, help="List of files to be used in the analysis, including raw repertoire files, metadata"
                                                   "files, additional files such as reference sequences or other inputs, existing dataset "
                                                   "files and others. Anything that immuneML can accept in the specification can be listed "
                                                   "here.")
    parser.add_argument("--tool", help="Name of the tool which calls immuneML. This name will be used to invoke appropriate API call, "
                                       "which will then preprocess the data/specs in tool-dependent way before running standard immuneML.")
    run_immuneML(parser.parse_args())


if __name__ == "__main__":
    main()
