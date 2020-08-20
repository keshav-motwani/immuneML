from collections import defaultdict

import numpy as np
import Levenshtein_search
from igraph import Graph


class NetworkConstructor:

    GRAPH_PROPERTIES = ['largest_component', "one_core", "max_k_core", "clique", "assortativity", "n_component", "transitivity", "density",
                        "diameter"]
    LOCAL_PROPERTIES = ["degree", "component_size", "transitivity", "authority", "pagerank", "eigenvector",
                        "betweenness", "closeness"]
    GRAPH_PROPERTIES = GRAPH_PROPERTIES + ["avg_" + i for i in LOCAL_PROPERTIES] + ["var_" + i for i in LOCAL_PROPERTIES]

    @staticmethod
    def construct_network(string_array: np.array, max_edit_distance: int = 1):

        string_to_index_map = NetworkConstructor.generate_string_to_index_map(string_array)

        strings = list(string_to_index_map.keys())
        wordset = Levenshtein_search.populate_wordset(-1, strings)

        result = []

        for index, string in enumerate(string_array):

            string_matches = Levenshtein_search.lookup(wordset, string, max_edit_distance)
            string_matches = [sequence[0] for sequence in string_matches]

            result.append(NetworkConstructor.get_edges(index, string_matches, string_to_index_map))

        result = [item for sublist in result for item in sublist]

        g = Graph(result).simplify()

        g.vs["sequence"] = string_array

        return g

    @staticmethod
    def generate_string_to_index_map(string_array):

        mapping = defaultdict(list)

        for index, string in enumerate(string_array):

            mapping[string].append(index)

        return mapping

    @staticmethod
    def get_edges(vertex_1_index, vertex_2_strings, string_to_index_map):

        result = [string_to_index_map[string] for string in vertex_2_strings]

        return [(vertex_1_index, item) for sublist in result for item in sublist]

    @staticmethod
    def compute_graph_properties(graph):

        N = graph.vcount()
        component_sizes = np.array(graph.components().sizes())
        degrees = np.array(graph.degree())

        local_transitivity = np.array(graph.transitivity_local_undirected())
        local_transitivity = local_transitivity[np.invert(np.isnan(local_transitivity))]
        authority = np.array(graph.authority_score())
        pagerank = np.array(graph.pagerank())
        eigenvector = np.array(graph.eigenvector_centrality())
        betweenness = 2 * np.array(graph.betweenness()) / ((N - 1) * (N - 2))
        closeness = np.array(graph.closeness(normalized=True))

        largest_component = 100 * max(component_sizes) / N
        one_core = 100 * np.sum(degrees > 0) / N
        max_k_core = 100 * degrees[degrees == np.max(degrees)].shape[0] / N
        clique = 100 * graph.omega() / N
        diameter = 100 * graph.diameter() / N
        assortativity = graph.assortativity_degree(directed=False)
        n_component = len(component_sizes)
        transitivity = graph.transitivity_undirected()
        density = 100 * graph.density()

        avg_degree = np.mean(graph.degree())
        avg_component_size = np.mean(component_sizes)
        avg_transitivity = np.mean(local_transitivity)
        avg_authority = np.mean(authority)
        avg_pagerank = np.mean(pagerank)
        avg_eigenvector = np.mean(eigenvector)
        avg_betweenness = np.mean(betweenness)
        avg_closeness = np.mean(closeness)

        var_degree = np.var(graph.degree())
        var_component_size = np.var(component_sizes)
        var_transitivity = np.var(local_transitivity)
        var_authority = np.var(authority)
        var_pagerank = np.var(pagerank)
        var_eigenvector = np.var(eigenvector)
        var_betweenness = np.var(betweenness)
        var_closeness = np.var(closeness)

        local_vars = locals()

        properties = [local_vars[i] for i in NetworkConstructor.GRAPH_PROPERTIES]

        data = {"property": NetworkConstructor.GRAPH_PROPERTIES,
                "value": np.array(properties)}

        return data
