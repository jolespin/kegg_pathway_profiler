import warnings
from collections import defaultdict, OrderedDict
from scipy.stats import false_discovery_control, hypergeom #, fisher_exact,nchypergeom_wallenius 
import pandas as pd
from tqdm import tqdm

from pyexeggutor import check_argument_choice
from .pathways import (
    update_graph_edge_weights_with_detected_kos,
    find_paths_in_pathway_graph,
)
# Feature Set Enrichment
def unweighted_set_enrichment(
    features:set, 
    feature_sets:dict, 
    background_set:set=None, 
    test_method:str="hypergeometric", 
    tol_test:float=0.05, 
    fdr_method:str="bh",
    ) -> pd.DataFrame: 
    """
    Perform unweighted set enrichment using Hypergeometric test (or Fisher Exact which is experimental)

    Args:
        features (set): Set of features to test
        feature_sets (dict): Feature sets with key as set set_name and value as feature set
        background_set (set, optional): Background feature set.  If None, then all features in `feature_sets` + `features` are used. [Defaults: None]
        test_method (str, optional): Type of test to use. Only {hypergeometric} is supported. [Defaults: "hypergeometric"]
        tol_test (float, optional): FDR threshold for signficance. [Defaults: 0.05]
        fdr_method (str, optional): FDR method 'bh' is for Benjamini-Hochberg, 'by' is for Benjaminini-Yekutieli. 
                                    The latter is more conservative, but it is guaranteed to control the FDR even 
                                    when the p-values are not from independent tests. [Defaults: "bh"]

    Raises:
        ValueError: features must be a subset of background_set

    Returns:
        _type_: _description_
        
        
    Future: 
     * Incorporate `feature_weights:pd.Series` using Wallenius' noncentral hypergeometric distribution

    Theory for hypergeometric test: 
    http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
    
    Source implementation: 
    https://github.com/jolespin/soothsayer
    
    """
    check_argument_choice(test_method, {"hypergeometric"})
    # if test_method == "fisher_exact":
    #     warnings.warn("test_method = 'fisher' is experimental")
    
    # Force set type for features in case other iterables are used
    features = set(features)
    
    # Force set type
    feature_sets = {k:set(v) for k,v in feature_sets.items()}

    # Union of all features in feature sets
    if not background_set:
        background_set = set.union(*feature_sets.values()) | features
    background_set = set(background_set)
    if not features <= background_set:
        raise ValueError(f"`features` must be a subset of `background_set` but `features` has {len(features - background_set)} unique features not in the background")
    
    number_of_features_in_background = len(background_set)
    
    data = defaultdict(OrderedDict)
    
    for set_name, set_features in feature_sets.items():
        # Get number of features in set
        number_of_features_in_set = len(set_features)
        # Get number of features in query
        number_of_features_in_query = len(features)
        # Get number of intersecting features between query and set
        intersecting_features = features & set_features
        number_of_features_intersection = len(intersecting_features)
        # Extra features
        extra_features = features - intersecting_features
        # Missing features
        # missing_features = set_features - intersecting_features
        
        # Rum hypergeometric test
        if test_method == "hypergeometric":
            model = hypergeom(
                M = number_of_features_in_background,
                n = number_of_features_in_set,
                N = number_of_features_in_query,
            )
            # "We want the *inclusive* upper tail : P-value = P(X≥x). 
            # For this, we can compute the exclusive upper tail of the value just below x. 
            # Indeed, since the distribution is discrete, P(X >x-1) = P(X ≥x)."
            # Source - http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
            p_value = model.sf(number_of_features_intersection - 1)
            data["p-value"][set_name] = p_value

        # if test_method == "fisher_exact": # The same as a hypergeometric test
        #     number_of_query_not_in_set = len(features - set_features)
        #     number_of_features_not_in_set = number_of_features_in_background - number_of_features_in_set
        #     contingency_table = [
        #         [number_of_features_intersection, number_of_features_in_query - number_of_features_intersection],
        #         [number_of_features_in_set - number_of_features_intersection, number_of_features_not_in_set - (number_of_features_in_query - number_of_features_intersection)],
        #     ]
        #     stat, p_value = fisher_exact(contingency_table, alternative="greater")
        #     data["Contingency Table"][set_name] = contingency_table
        #     data["Statistic"][set_name] = stat
        #     data["p-value"][set_name] = p_value

        # Store values
        data["Number Background Features (M)"][set_name] = number_of_features_in_background
        data["Number Set Features (n)"][set_name] = number_of_features_in_set
        data["Number Query Features (N)"][set_name] = number_of_features_in_query
        data["Number Intersecting Features (k)"][set_name] = number_of_features_intersection
        data["Intersecting Features"][set_name] = intersecting_features
        data["Extra Features"][set_name] = extra_features
        # data["Missing Features"][set_name] = missing_features


    # Create dataframe
    df = pd.DataFrame(data)
    df.insert(0, "Method", test_method)

    # Calculate adjusted p-value
    # df.insert(df.shape[1], "FDR", fdr_method)
    df.insert(df.shape[1], "FDR", false_discovery_control(df["p-value"], method=fdr_method))
    
    # Determine statistical significance
    if tol_test is not None:
        df.insert(df.shape[1], "FDR<{}".format(tol_test), df["FDR"] < tol_test)
        
    return df

def unweighted_pathway_enrichment_wrapper(
    evaluation_kos: set,
    database: dict,
    background_set: set = None,
    progressbar_description: str = "Calculating enrichment for most complete path in pathways:",
    **unweighted_set_enrichment_kws,
) -> pd.DataFrame:
    """
    Perform unweighted pathway enrichment analysis using the hypergeometric test on the most complete paths in pathways.

    This function identifies the most complete paths within pathways from a KEGG pathway database, based on the presence 
    of KEGG Orthology (KO) terms in the evaluation set. The function then performs a hypergeometric test to assess the 
    enrichment of these pathways, returning a DataFrame of enrichment results.

    Parameters
    ----------
    evaluation_kos : set
        A set of KEGG Orthology (KO) identifiers representing the features (e.g., genes) present in the dataset.
        
    database : dict
        A dictionary containing pathway data, where each key is a pathway identifier and the value is another dictionary 
        with pathway-specific data, including:
            - 'graph': A directed graph (as a networkx.MultiDiGraph) representing the pathway structure.
            - 'ko_to_nodes': A dictionary mapping KO identifiers to the nodes/edges in the graph that they represent.

    background_set : set, optional
        A set of KEGG Orthology (KO) identifiers representing the background set for the enrichment analysis. 
        If not provided, the union of all KO identifiers in the database will be used.

    progressbar_description : str, optional
        A description to display on the progress bar while processing pathways. Default is "Calculating enrichment for 
        most complete path in pathways:".

    **unweighted_set_enrichment_kws : dict
        Additional keyword arguments to pass to the `unweighted_set_enrichment` function, which performs the actual 
        hypergeometric test for enrichment.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the results of the pathway enrichment analysis. The DataFrame includes the 
        following columns:
            - 'p-value': The p-value from the hypergeometric test for enrichment.
            - 'odds_ratio': The odds ratio for enrichment.
            - 'intersecting_kos': The set of KO identifiers from the evaluation set that intersect with the pathway.
            - 'pathway_kos': The set of KO identifiers within the pathway.
            - Other columns may be added based on the output of the `unweighted_set_enrichment` function.

    Example
    -------
    >>> enrichment_results = unweighted_pathway_enrichment_wrapper(
    ...     evaluation_kos=my_kos, 
    ...     database=kegg_pathway_database,
    ...     background_set=all_kos
    ... )
    >>> print(enrichment_results.head())

    Notes
    -----
    This function first updates the edge weights in the pathway graphs based on the presence of KO identifiers in the 
    evaluation set. It then identifies the most complete paths within each pathway and calculates enrichment using the 
    hypergeometric test.
    """
    
    # Initialize a dictionary to store sets of KOs in the most complete paths for each pathway
    feature_sets = dict()
    
    # Iterate through each pathway in the database
    for id_pathway in tqdm(database, desc=progressbar_description, unit=" Pathways"):
        # Retrieve pathway data
        graph = database[id_pathway]["graph"]
        ko_to_nodes = database[id_pathway]["ko_to_nodes"]
    
        # Find intersecting KOs between the pathway and the evaluation set
        intersecting_kos = set(ko_to_nodes) & evaluation_kos
        
        if intersecting_kos:
            # Update the graph's edge weights based on detected KOs
            graph_weighted = update_graph_edge_weights_with_detected_kos(evaluation_kos, graph, ko_to_nodes)
            
            # Find the best path(s) based on the updated graph
            paths, path_to_ordered_kos, weights_normalized, most_complete_paths = find_paths_in_pathway_graph(graph_weighted)
        
            # Identify the most complete path and corresponding KOs
            most_complete_path = None
            kos_in_path = set()
            number_of_intersecting_ko = 0
            for path in most_complete_paths:
                query_kos = set(path_to_ordered_kos[path])
                n = len(query_kos & evaluation_kos)
                if n > number_of_intersecting_ko:
                    most_complete_path = path
                    kos_in_path = query_kos
                    number_of_intersecting_ko = n
            feature_sets[id_pathway] = kos_in_path

    # Perform unweighted set enrichment analysis on the most complete paths
    df_enrichment = unweighted_set_enrichment(
        features=evaluation_kos, 
        feature_sets=feature_sets, 
        background_set=background_set,
        **unweighted_set_enrichment_kws
    )
    
    # Set the index name for the resulting DataFrame
    df_enrichment.index.name = "id_pathway"
    
    return df_enrichment
