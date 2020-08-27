from collections import defaultdict

import numpy as np
import Levenshtein_search
from igraph import Graph


class NetworkConstructor:

    GRAPH_PROPERTIES = ['largest_component', "one_core", "max_k_core", "clique", "assortativity", "n_component", "transitivity", "density",]
                        # "diameter"]
    LOCAL_PROPERTIES = ["degree", "component_size", "transitivity", "authority", "pagerank", "eigenvector",]
                        # "betweenness", "closeness"]
    GRAPH_PROPERTIES = GRAPH_PROPERTIES + ["avg_" + i for i in LOCAL_PROPERTIES] + ["var_" + i for i in LOCAL_PROPERTIES]

    @staticmethod
    def construct_network(string_array: np.array, max_edit_distance: int = 1):

        strings = list(set(string_array))
        wordset = Levenshtein_search.populate_wordset(-1, strings)

        string_to_index_map = dict(zip(strings, range(len(strings))))

        result = []

        for index, string in enumerate(strings):

            string_matches = Levenshtein_search.lookup(wordset, string, max_edit_distance)
            index_matches = [string_to_index_map[sequence[0]] for sequence in string_matches]

            result.append([(index, match) for match in index_matches])

        result = [item for sublist in result for item in sublist]

        g = Graph(result).simplify()

        g.vs["sequence"] = strings

        return g

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
        # betweenness = 2 * np.array(graph.betweenness()) / ((N - 1) * (N - 2))
        # closeness = np.array(graph.closeness(normalized=True))

        largest_component = 100 * max(component_sizes) / N
        one_core = 100 * np.sum(degrees > 0) / N
        max_k_core = 100 * degrees[degrees == np.max(degrees)].shape[0] / N
        clique = 100 * graph.omega() / N
        # diameter = 100 * graph.diameter() / N
        assortativity = graph.assortativity_degree(directed=False)
        n_component = len(component_sizes)
        transitivity = graph.transitivity_undirected()
        density = 100 * graph.density()

        avg_degree = np.mean(degrees)
        avg_component_size = np.mean(component_sizes)
        avg_transitivity = np.mean(local_transitivity)
        avg_authority = np.mean(authority)
        avg_pagerank = np.mean(pagerank)
        avg_eigenvector = np.mean(eigenvector)
        # avg_betweenness = np.mean(betweenness)
        # avg_closeness = np.mean(closeness)

        var_degree = np.var(degrees)
        var_component_size = np.var(component_sizes)
        var_transitivity = np.var(local_transitivity)
        var_authority = np.var(authority)
        var_pagerank = np.var(pagerank)
        var_eigenvector = np.var(eigenvector)
        # var_betweenness = np.var(betweenness)
        # var_closeness = np.var(closeness)

        local_vars = locals()

        properties = [local_vars[i] for i in NetworkConstructor.GRAPH_PROPERTIES]

        data = {"property": NetworkConstructor.GRAPH_PROPERTIES,
                "value": np.array(properties)}

        return data
