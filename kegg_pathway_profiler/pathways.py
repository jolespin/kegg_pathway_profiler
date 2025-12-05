#!/usr/bin/env python
# from dataclasses import dataclass
import networkx as nx
from tqdm import tqdm
from pyexeggutor import (
    check_argument_choice,
    format_header,
)
from kegg_pathway_profiler.utils import (
    read_kos,

)
from kegg_pathway_profiler.parse import (
    parse_expression,
)


def update_graph_edge_weights_with_detected_kos(
    evaluation_kos: set, 
    graph: nx.MultiDiGraph, 
    ko_to_nodes: dict,
    ) -> nx.MultiDiGraph:
    """
    Updates the edge weights in the given directed graph based on the presence of KOs in the evaluation set.
    
    This function creates a copy of the provided graph to ensure the original is not modified. It then iterates 
    through each KO (KEGG Orthology identifier) in the evaluation set. For each KO, the function identifies the 
    corresponding edges in the graph and sets their 'weight_new' attribute to 0, indicating that these KOs have 
    been detected in the evaluation set.
    
    Parameters:
    ----------
    evaluation_kos : set
        A set of KEGG Orthology (KO) identifiers that have been detected or are under evaluation.
        
    graph : nx.MultiDiGraph
        A directed multigraph from the NetworkX library where nodes represent entities (e.g., genes, enzymes), 
        and edges represent interactions or relationships between these entities. Each edge is expected to have 
        a 'label' attribute corresponding to a KO and a 'weight_new' attribute that will be updated.
        
    ko_to_nodes : dict
        A dictionary mapping each KO identifier to a list of (start, end) node tuples, representing the edges in 
        the graph that are associated with that KO.
        
    Returns:
    -------
    nx.MultiDiGraph
        A modified copy of the input graph with updated edge weights based on the detected KOs.
    
    Notes:
    -----
    - The function assumes that each edge in the graph has a 'label' attribute matching a KO and a 'weight_new' 
      attribute that will be set to 0 for matching KOs.
    - If multiple edges exist between two nodes, the function finds and updates the first edge that matches the KO.
    
    Example:
    -------
    >>> evaluation_kos = {'K00100', 'K00200'}
    >>> ko_to_nodes = {'K00100': [(0, 1)], 'K00200': [(1, 2)]}
    >>> graph = nx.MultiDiGraph()
    >>> graph.add_edge(0, 1, label='K00100', weight=1, weight_new=1)
    >>> graph.add_edge(1, 2, label='K00200', weight=1, weight_new=1)
    >>> updated_graph = update_graph_edge_weights_with_detected_kos(evaluation_kos, graph, ko_to_nodes)
    >>> updated_graph[0][1][0]['weight_new']  # Should be 0
    >>> updated_graph[1][2][0]['weight_new']  # Should be 0
    """
    
    # Copy the graph to avoid modifying the original
    graph = graph.copy()
    
    # Update edge weights for the KOs present in the evaluation set
    for id_ko in evaluation_kos:
        if id_ko in ko_to_nodes:
            # Get the list of edges associated with the current KO
            edges = ko_to_nodes[id_ko]
            for start, finish in edges:
                if len(graph[start][finish]) > 0:
                    # Find the first edge that matches the KO label
                    for num in range(len(graph[start][finish])):
                        if graph[start][finish][num]['label'] == id_ko:
                            index = num
                            break
                else:
                    index = 0
                # Set the 'weight_new' attribute to 0 for edges that match the KO
                graph[start][finish][index]['weight_new'] = 0
                
    return graph

                
# Find pathways in pathway graph
def find_paths_in_pathway_graph(
    graph: nx.MultiDiGraph,
):
    """
    Finds and evaluates all paths in a directed, acyclic graph representing a biological pathway.

    This function identifies all possible paths through a graph where each node represents a biological
    component, and edges have associated weights. For each path, the function calculates a metric `M`, 
    defined as the ratio of new to old weights along the path. The function returns the paths, the sequence
    of labels associated with each path, the normalized weights (metric M), and the indices of the paths 
    with the smallest normalized weight.

    Parameters
    ----------
    graph : nx.MultiDiGraph
        A directed, acyclic graph where nodes represent biological components (e.g., genes, proteins) 
        and edges represent interactions or dependencies. Each edge should have associated attributes:
        'weight' (old weight) and 'weight_new' (new weight).

    Returns
    -------
    paths : list of lists
        A list where each sublist represents a sequence of nodes (as a path) in the graph.
        Example: [[0, 2, 4, 5, 6, 7, 1], [0, 3, 5, 7, 1]]
    
    path_to_ordered_kos : dict
        A dictionary where each key is a path index, and the value is a list representing the sequence 
        of labels (K numbers or identifiers) associated with the edges in that path.
        Example: {0: ['K00001', 'K00003', 'K00004', 'K00005'], 1: ['K00002', 'K00003', 'K00004', 'K00005']}
    
    weights_normalized : list of floats
        A list of normalized weights for each path, where each weight is the ratio of the new 
        to the old weight for that path.
        Example: [1.2, 1.1, 1.3]

    most_complete_paths : list of ints
        A list of indices of the paths that have the smallest normalized weight (metric M), 
        indicating the most "complete" or optimized paths.
        Example: [1, 2]

    Notes
    -----
    - The graph must be a directed acyclic graph (DAG) to ensure proper topological sorting.
    - The `graph` object is copied internally to avoid modifying the original graph.

    Example
    -------
    Suppose you have a pathway graph `G` where nodes are enzymes and edges represent reactions
    with weights representing reaction rates:

    >>> G = nx.MultiDiGraph()
    >>> # Add nodes and edges with weights
    >>> paths, path_to_ordered_kos, weights_normalized, most_complete_paths = find_paths_in_pathway_graph(G)
    >>> print(paths)
    [[0, 2, 4, 5, 6, 7, 1], [0, 3, 5, 7, 1]]
    >>> print(path_to_ordered_kos)
    {0: ['K00001', 'K00003', 'K00004', 'K00005'], 1: ['K00002', 'K00003', 'K00004', 'K00005']}
    >>> print(weights_normalized)
    [1.2, 1.1]
    >>> print(most_complete_paths)
    [1]
    """
    # graph = graph.copy()
    
    node_paths = {}
    path_labels = {}
    old_weights = {}
    new_weights = {}
    
    sorted_nodes = list(nx.topological_sort(graph))
    for node in sorted_nodes:
        record_count = 0
        node_paths[node] = []
        path_labels[node] = {}
        old_weights[node] = {}
        new_weights[node] = {}

        predecessors = graph.pred[node]
        if not predecessors:
            node_paths[node].append([node])
            path_labels[node][0] = []
            old_weights[node][0] = 0
            new_weights[node][0] = 0
            continue

        for pred in predecessors.keys():  # pred --> node
            number_of_paths_from_predecessor = len(path_labels[pred])

            for edge_key in predecessors[pred]:
                # 
                #    Handling multiple edges pred---A---->node
                #                        \____B____/ 
                for path_index in range(number_of_paths_from_predecessor):
                    current_labels = path_labels[pred][path_index]
                    path_labels[node][record_count] = \
                        current_labels + [predecessors[pred][edge_key]['label']]
                    old_weights[node][record_count] = \
                        old_weights[pred][path_index] + predecessors[pred][edge_key]['weight']
                    new_weights[node][record_count] = \
                        new_weights[pred][path_index] + predecessors[pred][edge_key]['weight_new']
                    record_count += 1
                for current_path in node_paths[pred]:
                    new_path = current_path + [node]
                    node_paths[node].append(new_path)

    paths = node_paths[1]
    path_to_ordered_kos = path_labels[1]
    old_weights_for_paths = old_weights[1]
    new_weights_for_paths = new_weights[1]

    weights_normalized = []
    for i in range(len(old_weights_for_paths)):
        weights_normalized.append(1.0 * new_weights_for_paths[i] / old_weights_for_paths[i])

    min_weight_normalized = min(weights_normalized)  # Find the minimum metric value
    most_complete_paths = []  # Initialize an empty list to store indices

    # Iterate through the metrics list to find indices with the minimum metric
    for i, w in enumerate(weights_normalized):
        if w == min_weight_normalized:
            most_complete_paths.append(i)

    return paths, path_to_ordered_kos, weights_normalized, most_complete_paths


def get_pathway_coverage(
    evaluation_kos: set, 
    graph: nx.MultiDiGraph, 
    ko_to_nodes: dict, 
    optional_kos: set,
):
    """
    Calculate the coverage of a pathway graph based on a given set of KOs (KEGG Orthology terms).

    This function evaluates the extent to which a set of KOs is represented in a pathway graph. It modifies 
    the graph by updating edge weights based on the presence of KOs in the `evaluation_kos` set. The coverage 
    is determined by finding the most complete path in the graph and calculating the proportion of KOs that 
    are covered.

    Parameters
    ----------
    evaluation_kos : set
        A set of KOs to evaluate against the pathway graph.
    
    graph : nx.MultiDiGraph
        A directed, acyclic graph representing the pathway. Nodes represent biological components, and edges 
        represent interactions or dependencies, with weights indicating the strength or presence of these 
        connections.
    
    ko_to_nodes : dict
        A dictionary mapping each KO to its corresponding edges in the graph. Each KO is associated with a list 
        of tuples representing the start and end nodes of the edges it influences.
    
    optional_kos : set
        A set of KOs that are optional and not required for full pathway coverage.

    Returns
    -------
    dict
        A dictionary containing the following keys:
        - 'coverage': float, the proportion of the pathway covered by the given KOs (0 to 1 scale).
        - 'number_of_best_paths': int, the number of paths with the highest coverage.
        - 'required_kos_in_path': set, KOs from `evaluation_kos` that are present in the most complete path.
        - 'required_kos_missing_in_path': set, KOs from `evaluation_kos` that are missing from the most complete path.

    Notes
    -----
    - The function modifies a copy of the input graph, leaving the original graph unchanged.
    - The function assumes that the graph is directed and acyclic, suitable for pathway analysis.
    """
    # Updates the edge weights in the given directed graph based on the presence of KOs in the evaluation set.
    graph_weighted = update_graph_edge_weights_with_detected_kos(evaluation_kos, graph, ko_to_nodes)
    
    # Find the best path(s) based on the updated graph
    paths, path_to_ordered_kos, weights_normalized, most_complete_paths = find_paths_in_pathway_graph(graph_weighted)

    # Select the first of the most complete paths (all have the same coverage)
    most_complete_path = most_complete_paths[0]
    coverage = 1 - weights_normalized[most_complete_path]

    # Determine which KOs are required but missing in the most complete path
    required_kos_missing_in_path = set()
    required_kos_in_path = set()
    if coverage > 0:
        kos_in_path = set(path_to_ordered_kos[most_complete_path])
        kos_missing_in_path = evaluation_kos - kos_in_path
        required_kos_missing_in_path = kos_missing_in_path - optional_kos
        required_kos_in_path = kos_in_path & evaluation_kos

    # Return a dictionary containing coverage information
    return dict(
        coverage=coverage,
        number_of_best_paths=len(most_complete_paths),
        most_complete_path=path_to_ordered_kos[most_complete_path],
        required_kos_in_path=required_kos_in_path,
        required_kos_missing_in_path=required_kos_missing_in_path,
    )

def pathway_coverage_wrapper(
    evaluation_kos: set,
    database: dict,
    progressbar_description: str = "Calculating pathway coverage:",
    progressbar = False,
) -> dict:
    """
    Calculates the coverage of pathways in a KEGG database based on a set of evaluation KOs (KEGG Orthology identifiers).
    
    This function iterates over pathways in a given KEGG database, computes the coverage of each pathway by 
    comparing its KOs with the provided evaluation set, and stores the results in a dictionary. The progress 
    of the computation is displayed using a progress bar.
    
    Parameters:
    ----------
    evaluation_kos : set
        A set of KEGG Orthology (KO) identifiers that are to be evaluated for pathway coverage.
        
    database : dict
        A dictionary representing the KEGG database, where each key is a pathway identifier and the value is 
        a dictionary containing pathway data. This data should include:
        - "graph": A NetworkX directed multigraph representing the pathway.
        - "ko_to_nodes": A dictionary mapping each KO to a list of (start, end) node tuples in the pathway graph.
        - "optional_kos": A set of optional KOs that are not required for pathway coverage.
    
    progressbar_description : str, optional
        A string to describe the progress bar shown during the calculation, by default "Calculating pathway coverage:".
        
    Returns:
    -------
    dict
        A dictionary where each key is a pathway identifier, and the value is the result of the `get_pathway_coverage` 
        function, which includes coverage statistics and other relevant information.
    
    Example:
    -------
    >>> evaluation_kos = {'K00100', 'K00200', 'K00300'}
    >>> database = {
    ...     'pathway1': {
    ...         'graph': nx.MultiDiGraph(),
    ...         'ko_to_nodes': {'K00100': [(0, 1)], 'K00200': [(1, 2)]},
    ...         'optional_kos': {'K00300'}
    ...     },
    ...     'pathway2': {
    ...         'graph': nx.MultiDiGraph(),
    ...         'ko_to_nodes': {'K00400': [(2, 3)]},
    ...         'optional_kos': set()
    ...     }
    ... }
    >>> results = pathway_coverage_wrapper(evaluation_kos, database)
    >>> results  # Will contain coverage information for 'pathway1' and 'pathway2'
    """
    
    pathway_to_results = dict()  # Dictionary to store coverage results for each pathway
    
    # Iterate over each pathway in the database
    if progressbar:
        iterable = tqdm(database, desc=progressbar_description, unit=" Pathways")
    else:
        iterable = database
    for id_pathway in iterable:
        # Extract the graph, KO-to-nodes mapping, and optional KOs for the current pathway
        graph = database[id_pathway]["graph"]
        ko_to_nodes = database[id_pathway]["ko_to_nodes"]
        optional_kos = database[id_pathway]["optional_kos"]
        
        # Find the intersection of evaluation KOs and the KOs present in the pathway
        intersecting_kos = set(ko_to_nodes) & evaluation_kos
        
        # If there are intersecting KOs, calculate pathway coverage
        if intersecting_kos:
            results = get_pathway_coverage(
                evaluation_kos=evaluation_kos,
                graph=graph, 
                ko_to_nodes=ko_to_nodes, 
                optional_kos=optional_kos,
            )
            pathway_to_results[id_pathway] = results  # Store the results in the dictionary
    
    return pathway_to_results


# @dataclass
class Pathway:
    """
    A class for processing KEGG pathway definitions to generate a MultiDiGraph and associated mappings.
    
    Attributes:
        id (str): The identifier for the pathway.
        definition (str): The KEGG pathway definition string.
        name (str, optional): The name of the pathway. Defaults to None.
        classes (str, optional): The classes associated with the pathway. Defaults to None.
        start_node_color (str): Color of the start node in the graph. Defaults to "green".
        end_node_color (str): Color of the end node in the graph. Defaults to "red".
    """
    def __init__(
        self,
        id: str,
        definition: str,
        name: str = None,
        classes: str = None,
        start_node_color: str = "green",
        end_node_color: str = "red",
        ):
        self.id = id
        self.definition = definition 
        self.name = name 
        self.classes = classes 
        self.start_node_color = start_node_color 
        self.end_node_color = end_node_color
        
        # Build graph
        graph = nx.MultiDiGraph(
            id=self.id, 
            definition=self.definition,
            name=self.name, 
            classes=self.classes, 
        )
        graph.add_node(0, color=self.start_node_color)
        graph.add_node(1, color=self.end_node_color)

        # Parse definition
        self.graph_, self.ko_to_nodes_, self.optional_kos_ = parse_expression(
            graph=graph,
            ko_to_nodes={},
            optional_kos=set(),
            expression=self.definition,
            start_node=0,
            end_node=1,
            weight=1
        )
        self.kos_ = set(self.ko_to_nodes_)
        
    # Evaluate pathway coverage from set of KOs

    def evaluate(
        self,
        kos:set, 
        return_type:str=None,
        ):
        """  Evaluate pathway coverage from set of KOs


        Args:
            kos (set): _description_
            return_type (str, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        choices = {
            None, "coverage", "number_of_best_paths", "most_complete_path", "required_kos_in_path", "required_kos_missing_in_path",
        }
        check_argument_choice(return_type, choices)
        results = get_pathway_coverage(
            evaluation_kos=kos,
            graph=self.graph_, 
            ko_to_nodes=self.ko_to_nodes_, 
            optional_kos=self.optional_kos_,
        )

        if return_type:
            return results[return_type]
        else:
            return results
                        
    def __repr__(self):
        name_text = f"{self.__class__.__name__}(id:{self.id})"
        n = len(name_text)
        pad = 4
        fields = [
            format_header(name_text,line_character="=", n=n),        
            "Properties:",
            pad*" " + f"- name: {self.name}",
            pad*" " + f"- classes: {self.classes}",
            pad*" " + f"- number_of_kos: {len(self.kos_)}",
            "Definition:",
            pad*" " + self.definition,
            ]

        return "\n".join(fields)


def pathway_graph_wrapper(
    pathway_to_definition, 
    description="Parsing pathway definitions", 
    return_type=tuple
    ):
    """
    Parse pathway definitions and create Pathway objects or graph data.

    This function iterates over a dictionary of pathway definitions, creates a `Pathway` 
    object for each definition, and stores the processed data. The return type can be 
    either a `Pathway` object or a tuple containing graph-related data.

    Args:
        pathway_to_definition (dict): A dictionary where each key is a pathway ID and the value 
                                      is the corresponding KEGG pathway definition string.
        description (str, optional): A description for the progress bar. Defaults to "Parsing pathway definitions".
        return_type (type, optional): The type of object to return for each pathway. Can be `tuple` or `Pathway`.
                                      Defaults to `tuple`.

    Returns:
        dict: A dictionary mapping each pathway ID to either a `Pathway` object or a tuple 
              containing the graph, KO-to-node mapping, and optional KOs set.
    """
    check_argument_choice(return_type, {tuple, Pathway})
    
    pathway_to_data = {}
    for id_pathway, definition in tqdm(pathway_to_definition.items(), desc=description, unit=" Pathways"):
        pathway = Pathway(id=id_pathway, definition=definition)
        pathway_to_data[id_pathway] = pathway

    if return_type == tuple:
        for id_pathway, pathway in pathway_to_data.items():
            data = (pathway.graph_, pathway.ko_to_nodes_, pathway.optional_kos_)
            pathway_to_data[id_pathway] = data

    return pathway_to_data
        
        
