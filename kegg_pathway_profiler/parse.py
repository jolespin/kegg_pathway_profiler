#!/usr/bin/env python
import networkx as nx


def order_separators(level_to_positions: dict) -> list:
    """
    Determine the correct parsing order for separators based on their levels.

    The separators are sorted and processed in the order of [',', ' ', '+', '-'].

    Args:
        level_to_positions (dict): A dictionary mapping separator levels to their positions in the expression.

    Returns:
        list: A list of separators in the correct parsing order.
    """
    if level_to_positions:
        # Extract levels from the keys and determine the range of levels
        levels = sorted(level_to_positions)
        min_level = int(levels[0].split('_')[0])
        max_level = int(levels[-1].split('_')[0]) + 1

        # Create the ordered separators based on levels and predefined order
        ordered_separators = [f"{level}_{separator}" for level in range(min_level, max_level)
                              for separator in [",", " ", "+", "-"]]

        # Filter and return separators present in the provided levels
        return [separator for separator in ordered_separators if separator in levels]
    
    return []


def update_level_to_positions(level_to_positions: dict, symbol: str, level: int, index: int) -> dict:
    """
    Update the dictionary with the positions of spaces or commas based on their level.

    Args:
        level_to_positions (dict): The current dictionary mapping levels to positions.
        symbol (str): The symbol (space or comma) to update positions for.
        level (int): The current level in the expression.
        index (int): The position of the symbol in the expression.

    Returns:
        dict: Updated dictionary with the new positions added.
    """
    key = f"{level}_{symbol}"
    
    # Initialize the list for the key if it doesn't exist
    if key not in level_to_positions:
        level_to_positions[key] = []
    
    # Append the index to the corresponding level and symbol
    level_to_positions[key].append(index)
    return level_to_positions


def get_bracket_levels(expression: str) -> list:
    """
    Determine the nesting levels of brackets in the expression.

    Args:
        expression (str): The input expression containing brackets.

    Returns:
        list: A list indicating the level of each character in the expression.
    """
    levels = []
    current_level_stack = []
    current_level = -1

    # Iterate through each character to calculate bracket levels
    for character in expression:
        if character == '(':
            current_level += 1
            current_level_stack.append(current_level)
            levels.append(current_level)
        elif character == ')':
            levels.append(current_level_stack.pop())
        else:
            levels.append(-1)
    
    return levels


def parse_levels(expression: str) -> dict:
    """
    Create a dictionary of separators and their positions within the expression.

    The keys are formatted as 'level_separator' (e.g., '1_,' or '0_ '), and the values
    are lists of positions where these separators occur.

    Args:
        expression (str): The expression to parse.

    Returns:
        dict: A dictionary mapping separators to their positions within the expression.
    """
    level_to_positions = {}
    current_level = 0
    index = 0

    # Iterate through the expression to identify and record separator positions
    while index < len(expression):
        character = expression[index]
        if character in {' ', ',', '-', '+'}:
            level_to_positions = update_level_to_positions(level_to_positions, character, current_level, index)
        elif character == '(':
            current_level += 1
        elif character == ')':
            current_level -= 1
        else:
            # Skip over continuous non-separator characters
            while index < len(expression) - 1 and expression[index + 1] not in {' ', ',', '(', ')', '-', '+'}:
                index += 1
        index += 1
    
    return level_to_positions


def strip_outer_brackets(expression: str, bracket_levels: list) -> str:
    """
    Remove outer brackets if the expression is wrapped in them.

    Args:
        expression (str): The expression to strip brackets from.
        bracket_levels (list): The levels of brackets in the expression.

    Returns:
        str: The expression with outer brackets removed, if applicable.
    """
    if expression.startswith('(') and expression.endswith(')') and bracket_levels[0] == bracket_levels[-1]:
        return expression[1:-1]
    return expression


def parse_expression(graph: nx.MultiDiGraph, ko_to_nodes: dict, optional_kos: set, expression: str,
                     start_node: int, end_node: int, weight: float) -> tuple:
    """
    Recursively parse the expression to build the graph structure and update KO mappings.

    Args:
        graph (nx.MultiDiGraph): The graph being constructed.
        ko_to_nodes (dict): Mapping of KO identifiers to the corresponding nodes in the graph.
        optional_kos (set): Set of optional KO identifiers.
        expression (str): The expression to parse.
        start_node (int): The start node in the graph.
        end_node (int): The end node in the graph.
        weight (float): The weight of the edge being added.

    Returns:
        tuple: Updated graph, KO-to-node mapping, and optional KOs set.
    """
    # Handle the case of missing KO with a placeholder
    if expression == '--':
        missing_ko = 'K00000'
        graph.add_edge(start_node, end_node, label=missing_ko, weight=0, weight_new=0, name='-')
        optional_kos.add(missing_ko)
        ko_to_nodes.setdefault(missing_ko, []).append([start_node, end_node])
        return graph, ko_to_nodes, optional_kos

    # Strip outer brackets and parse the expression by levels
    expression = strip_outer_brackets(expression, get_bracket_levels(expression))
    level_to_positions = parse_levels(expression)
    separator_order = order_separators(level_to_positions)

    # Handle single optional KO case
    if len(separator_order) == 1 and separator_order[0] == '0_-' and expression.startswith('-'):
        ko_id = expression[1:]
        graph.add_edge(start_node, end_node, label=ko_id, weight=0, weight_new=0, name='-')
        optional_kos.add(ko_id)
        ko_to_nodes.setdefault(ko_id, []).append([start_node, end_node])
        return graph, ko_to_nodes, optional_kos

    # Recursively parse sub-expressions based on separator order
    if separator_order:
        separator_key = separator_order[0]
        separator_symbol = separator_key.split('_')[1]
        sub_weight = weight / (len(level_to_positions[separator_key]) + 1) if separator_symbol in {'+', ' '} else weight

        current_separator = 0
        current_start_node = start_node
        current_end_node = end_node

        # Process each subexpression divided by the separator
        for separator_index in sorted(level_to_positions[separator_key]):
            if separator_symbol in {' ', '+', '-'}:
                current_end_node = len(graph)
                graph.add_node(current_end_node)

            subexpression = expression[current_separator:separator_index]
            graph, ko_to_nodes, optional_kos = parse_expression(
                graph, ko_to_nodes, optional_kos, subexpression, current_start_node, current_end_node, sub_weight
            )
            current_separator = separator_index + 1

            if separator_symbol in {' ', '+', '-'}:
                current_start_node = current_end_node

        # Process the remaining expression after the last separator
        remaining_expression = expression[current_separator:]
        if separator_symbol in {' ', '+', '-'}:
            current_start_node = current_end_node
            current_end_node = end_node

        if separator_symbol == '-' and level_to_positions[separator_key]:
            sub_weight = 0

        graph, ko_to_nodes, optional_kos = parse_expression(
            graph, ko_to_nodes, optional_kos, remaining_expression, current_start_node, current_end_node, sub_weight
        )
        return graph, ko_to_nodes, optional_kos
    else:
        # Base case: add a single edge for the terminal expression
        graph.add_edge(start_node, end_node, label=expression, weight=weight, weight_new=weight, name='node')
        ko_to_nodes.setdefault(expression, []).append([start_node, end_node])
        if weight == 0:
            optional_kos.add(expression)
        return graph, ko_to_nodes, optional_kos


# def process_pathway(id_pathway: str, definition: str, start_node_color: str = "green",
#                     end_node_color: str = "red") -> tuple:
#     """
#     Process a single pathway definition to generate a graph and KO-to-node mappings.

#     Args:
#         id_pathway (str): The identifier for the pathway being processed.
#         definition (str): The KEGG pathway definition string.
#         start_node_color (str, optional): Color of the start node. Defaults to "green".
#         end_node_color (str, optional): Color of the end node. Defaults to "red".

#     Returns:
#         tuple: A tuple containing the graph, KO-to-node mapping, and optional KOs set.
#     """
#     graph = nx.MultiDiGraph(name=id_pathway)
    
#     # Add start and end nodes with specified colors
#     graph.add_node(0, color=start_node_color)
#     graph.add_node(1, color=end_node_color)

#     graph, ko_to_nodes, optional_kos = parse_expression(
#         graph=graph,
#         ko_to_nodes={},
#         optional_kos=set(),
#         expression=definition,
#         start_node=0,
#         end_node=1,
#         weight=1
#     )
#     return graph, ko_to_nodes, optional_kos

