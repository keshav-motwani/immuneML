from unittest import TestCase

from source.dsl.MLParser import MLParser
from source.dsl.SymbolTable import SymbolTable
from source.dsl.SymbolType import SymbolType
from source.ml_methods.LogisticRegression import LogisticRegression


class TestMLParser(TestCase):
    def test_parse_ml_methods(self):

        params = {
            "ml_methods": {
                "LR1": {
                    "type": "LogisticRegression",
                    "params": {
                        "max_iter": 1000,
                        "penalty": "l1",
                    },
                    "encoding": "e1",
                    "labels": ["CD"],
                    "metrics": ["accuracy", "balanced_accuracy"],
                    "min_example_count": 1
                },
                "LR2": {
                    "type": "LogisticRegression",
                    "encoding": "e1",
                    "labels": ["CD"],
                    "metrics": ["accuracy", "balanced_accuracy"],
                    "min_example_count": 1
                },
                "SVM1": {
                    "type": "SVM",
                    "params": {
                        "max_iter": [1000, 2000],
                        "penalty": ["l1", "l2"]
                    },
                    "encoding": "e1",
                    "labels": ["CD"],
                    "metrics": ["accuracy", "balanced_accuracy"],
                    "split_count": 1,
                    "model_selection_cv": False,
                    "model_selection_n_folds": -1,
                    "assessment_type": "LOOCV",
                }
            }
        }

        symbol_table = SymbolTable()
        with self.assertRaises(AssertionError):
            MLParser.parse(params, symbol_table)
        symbol_table.add("e1", SymbolType.ENCODING, {})
        symbol_table, desc = MLParser.parse(params, symbol_table)
        self.assertTrue(symbol_table.get("SVM1")["method"]._parameter_grid is not None and len(symbol_table.get("SVM1")["method"]._parameter_grid["max_iter"]) == 2)
        self.assertTrue(symbol_table.get("LR1")["method"]._parameters is not None and symbol_table.get("LR1")["method"]._parameters["penalty"] == "l1")
        self.assertTrue(isinstance(symbol_table.get("LR2")["method"], LogisticRegression))

        self.assertEqual("SVM", desc["SVM1"]["type"])
        self.assertEqual(2, desc["SVM1"]["min_example_count"])
        self.assertEqual("random", desc["LR1"]["assessment_type"])
        self.assertEqual(5, desc["LR1"]["split_count"])