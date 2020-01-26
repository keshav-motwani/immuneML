from unittest import TestCase
import yaml
import os
import numpy as np
import pandas as pd

from source.environment.EnvironmentSettings import EnvironmentSettings
from source.ml_methods.SimpleLogisticRegression import SimpleLogisticRegression
from source.reports.ml_reports.LogisticRegressionCoefficients import LogisticRegressionCoefficients
from source.util.PathBuilder import PathBuilder

class TestLogisticRegressionCoefficients(TestCase):

    def _create_dummy_lr_model(self, path):
        # dummy logistic regression with 100 observations with 20 features belonging to 2 classes
        dummy_lr = SimpleLogisticRegression()
        dummy_lr.fit_by_cross_validation(np.random.rand(100, 20),
                                         {"l1": [i % 2 for i in range(0, 100)]},
                                         number_of_splits=2,
                                         label_names=["l1"])

        # Change coefficients to values 1-20
        dummy_lr.models["l1"].coef_ = np.array(list(range(0, 20))).reshape(1, -1)

        with open(path + "ml_details.yaml", "w") as file:
            yaml.dump({"l1": {"feature_names": [f"feature{i}" for i in range(20)]}},
                      file)

        return dummy_lr

    def _create_report(self, path):
        report = LogisticRegressionCoefficients(**{"coefs_to_plot": ["all", "nonzero", "cutoff", "n_largest"],
                                                   "cutoff": [10],
                                                   "n_largest": [5]})

        report.method = self._create_dummy_lr_model(path)
        report.ml_details_path = path + "ml_details.yaml"
        report.label = "l1"
        report.result_path = path

        return report

    def test_generate(self):
        path = EnvironmentSettings.root_path + "test/tmp/logregcoefsreport/"
        PathBuilder.build(path)

        report = self._create_report(path)

        # Running the report
        report.check_prerequisites()
        report.generate()

        # Actual tests
        self.assertTrue(os.path.isfile(path + "coefficients.csv"))
        self.assertTrue(os.path.isfile(path + "all_coefficients.pdf"))
        self.assertTrue(os.path.isfile(path + "nonzero_coefficients.pdf"))
        self.assertTrue(os.path.isfile(path + "cutoff_10_coefficients.pdf"))
        self.assertTrue(os.path.isfile(path + "largest_5_coefficients.pdf"))

        written_data = pd.read_csv(path + "coefficients.csv")

        self.assertListEqual(list(written_data.columns), ["coefficients", "features"])
        self.assertListEqual(list(written_data["coefficients"]), [i for i in range(20)])
        self.assertListEqual(list(written_data["features"]), [f"feature{i}" for i in range(20)])

