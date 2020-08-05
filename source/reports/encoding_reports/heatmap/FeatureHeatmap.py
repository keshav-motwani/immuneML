import copy
import json

import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP

pandas2ri.activate()

from source.data_model.dataset.RepertoireDataset import RepertoireDataset
from source.reports.encoding_reports.EncodingReport import EncodingReport
from source.environment.EnvironmentSettings import EnvironmentSettings
from source.util.PathBuilder import PathBuilder


class FeatureHeatmap(EncodingReport):

    """
    Generate annotated heatmap of encoded data matrix
    @param feature_annotations: list of column names to annotate from encoded_dataset.encoded_data.feature_annotations -
    this varies based on the encoding performed, so these column names should be chosen correctly
    @param example_annotations: list of keys to annotate from encoded_dataset.encoded_data.labels - these are added to
    from the metadata file based on what labels were encoded (based on LabelConfiguration)
    @param one_hot_encode_feature_annotations: same as feature_annotations, but only discrete columns - these columns
    are then one-hot encoded and annotated
    @param one_hot_encode_example_annotations: same as one_hot_encode_feature_annotations, but for example annotations
    @param palette: dictionary with each annotation (from both feature_annotations and example_annotations) that you
    want to specify particular color palettes for - it should be in the following format:
        {
            "annotation_1": {"A": "blue", "B": "red", ...} if discrete
            "annotation_2": {"colors": ["blue", "white", "red"], "breaks": [-1, 0, 1]} if continuous
        }
    if it's discrete, and all values do not have a color specified, default values will be filled in, and if it's a
    continuous column, if only colors are defined, then breaks will be equally spaced between the min and max values
    for that column
    @param cluster_features: boolean whether to cluster features (rows)
    @param cluster_examples: boolean whether to cluster examples (column)
    @param subset_nonzero_features: boolean whether to remove features with only 0 values
    @param show_feature_dend: boolean whether to show the dendrogram for hierarchical clustering of features
    @param show_example_dend: boolean whether to show the dendrogram for hierarchical clustering of examples
    @param show_feature_names: boolean whether to show feature names
    @param show_example_names: boolean whether to show example names
    @param show_legend_features: if all legends are not desired, then a subset of feature_annotations can be specified
    here to specify which legends for feature annotations are desire
    @param show_legend_examples: same as show_legend_features but for example_annotations
    @param legend_position: either "bottom" or "side" specifying where all legends should be drawn
    @param text_size: text size for global text on heatmap
    @param feature_names_size: text size for feature names
    @param example_names_size: text size for example names
    @param scale_features: boolean whether to min-max scale features
    @param height: float indicating final figure height
    @param width: float indicating final figure width
    @param result_name: string indicating resulting figure file name
    @param result_path: string indicating resulting figure file path

    example:

    one_hot_encode_example_annotations=["disease"],
    example_annotations=["age", "week"],
    feature_annotations=["antigen"],
    palette={"week": {"0": "#BE9764"}, "antigen": {"GAD": "cornflowerblue", "INSB": "firebrick"},
             "age": {"colors": ["blue", "white", "red"], "breaks": [0, 20, 100]}},
    result_path=path,
    show_feature_names=True,
    feature_names_size=7,
    show_example_names=True,
    example_names_size=1,
    scale_features=False,
    text_size=9,
    height=6,
    width=6
    """

    FEATURE = "feature"
    EXAMPLE = "example"

    @classmethod
    def build_object(cls, **kwargs):
        return FeatureHeatmap(**kwargs)

    def __init__(self, dataset: RepertoireDataset = None, feature_annotations: list = [],
                 example_annotations: list = [], one_hot_feature_annotations: list = [],
                 one_hot_example_annotations: list = [], palette: dict = {},
                 cluster_features: bool = True, cluster_examples: bool = True, show_feature_dend: bool = True,
                 show_example_dend: bool = True, show_feature_names: bool = False, show_example_names: bool = False,
                 show_legend_features: list = None, show_legend_examples: list = None,
                 legend_position: str = "side", text_size: float = 10, feature_names_size: float = 7, example_names_size: float = 7,
                 feature_standardization: str = "NULL", lower_quantile: float = 0.01, upper_quantile: float = 0.99,
                 heatmap_color: str = "BWR", color_by_quantile=True, height: float = 10,
                 width: float = 10, result_name: str = "feature_heatmap",
                 result_path: str = None):

        super().__init__()
        self.dataset = dataset
        self.feature_annotations = feature_annotations
        self.example_annotations = example_annotations
        self.one_hot_feature_annotations = one_hot_feature_annotations
        self.one_hot_example_annotations = one_hot_example_annotations
        self.palette = palette
        self.cluster_features = cluster_features
        self.cluster_examples = cluster_examples
        self.show_feature_dend = show_feature_dend
        self.show_example_dend = show_example_dend
        self.show_feature_names = show_feature_names
        self.show_example_names = show_example_names
        self.show_legend_features = show_legend_features
        self.show_legend_examples = show_legend_examples
        self.legend_position = legend_position
        self.text_size = text_size
        self.feature_names_size = feature_names_size
        self.example_names_size = example_names_size
        self.feature_standardization = feature_standardization
        self.lower_quantile = lower_quantile
        self.upper_quantile = upper_quantile
        self.heatmap_color = heatmap_color
        self.color_by_quantile = color_by_quantile
        self.height = height
        self.width = width
        self.result_name = result_name
        self.result_path = result_path

        if self.show_legend_features is None:
            self.show_legend_features = copy.deepcopy(self.feature_annotations)

        if self.show_legend_examples is None:
            self.show_legend_examples = copy.deepcopy(self.example_annotations)

    def generate(self):
        PathBuilder.build(self.result_path)
        self._safe_plot(output_written=False)

    def _plot(self):

        matrix = self._prepare_matrix()
        feature_annotations = self._prepare_annotations(self.dataset.encoded_data.feature_annotations,
                                                        FeatureHeatmap.FEATURE)
        example_annotations = self._prepare_annotations(pd.DataFrame(self.dataset.encoded_data.labels),
                                                        FeatureHeatmap.EXAMPLE)
        one_hot_feature_annotations = self._prepare_one_hot_annotations(self.dataset.encoded_data.feature_annotations,
                                                                        FeatureHeatmap.FEATURE)
        one_hot_example_annotations = self._prepare_one_hot_annotations(pd.DataFrame(self.dataset.encoded_data.labels),
                                                                        FeatureHeatmap.EXAMPLE)

        self._save_results(matrix, feature_annotations, example_annotations, one_hot_feature_annotations, one_hot_example_annotations)

        with open(EnvironmentSettings.root_path + "source/visualization/Heatmap.R") as f:
            string = f.read()

        plot = STAP(string, "plot")

        plot.plot_heatmap(matrix=matrix,
                          row_annotations=feature_annotations,
                          column_annotations=example_annotations,
                          one_hot_row_annotations=one_hot_feature_annotations,
                          one_hot_column_annotations=one_hot_example_annotations,
                          palette=json.dumps(self.palette),
                          row_names=list(self.dataset.encoded_data.feature_names),
                          column_names=list(self.dataset.encoded_data.example_ids),
                          cluster_rows=self.cluster_features,
                          cluster_columns=self.cluster_examples,
                          show_row_dend=self.show_feature_dend,
                          show_column_dend=self.show_example_dend,
                          show_row_names=self.show_feature_names,
                          show_column_names=self.show_example_names,
                          show_legend_row=self.show_legend_features,
                          show_legend_column=self.show_legend_examples,
                          legend_position = self.legend_position,
                          text_size=self.text_size,
                          row_names_size=self.feature_names_size,
                          column_names_size=self.example_names_size,
                          row_standardization=self.feature_standardization,
                          lower_quantile=self.lower_quantile,
                          upper_quantile=self.upper_quantile,
                          heatmap_color=self.heatmap_color,
                          color_by_quantile=self.color_by_quantile,
                          height=self.height,
                          width=self.width,
                          result_path=self.result_path,
                          result_name=self.result_name)

    def _prepare_matrix(self):
        matrix = self.dataset.encoded_data.examples.A.T
        return matrix

    def _prepare_annotations(self, data, type):
        return data[getattr(self, type + "_annotations")]

    def _prepare_one_hot_annotations(self, data, type):
        return data[getattr(self, "one_hot_" + type + "_annotations")]

    def _save_results(self, matrix):

        np.savetxt(self.result_path + "/" + self.result_name + "_matrix.csv", matrix, delimiter=",")

        self.dataset.encoded_data.feature_annotations.to_csv(self.result_path + "/" + self.result_name + "_feature_annotations.csv")

        pd.DataFrame(self.dataset.encoded_data.labels).to_csv(self.result_path + "/" + self.result_name + "_example_annotations.csv")
