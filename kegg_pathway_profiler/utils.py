#!/usr/bin/env python
from typing import Union, Set
from collections import defaultdict
from tqdm import tqdm
from pyexeggutor import open_file_reader

def read_kos(filepath: str, name: str):
    """
    Reads KEGG Orthology (KO) identifiers from a file and maps them to corresponding genomes.

    This function processes a file containing KO identifiers, either for a single genome or multiple genomes.
    The file can either have one KO per line or be a tab-separated table with two columns: 
    [id_genome] <tab> [id_ko]. It returns a dictionary mapping genome identifiers to sets of KO identifiers.

    Parameters
    ----------
    filepath : str
        The path to the input file containing KO identifiers.
    
    name : str
        The identifier for the genome if the file contains KOs for a single genome.

    Returns
    -------
    genome_to_kos : dict
        A dictionary where keys are genome identifiers and values are sets of KO identifiers.
    
    Raises
    ------
    IndexError
        If the input file has a number of columns that is not 1 or 2.
    
    ValueError
        If any KO identifier in the file is not 6 characters long or does not start with 'K'.
    
    Notes
    -----
    - The function assumes that the input file has no header.
    - The function checks the format of KO identifiers and raises an error if they do not conform to the expected format.
    """
    
    # Initialize a dictionary to store genome-to-KO mappings
    genome_to_kos = defaultdict(set)
    
    with open_file_reader(filepath) as f:
        # Read the first line and determine the number of columns
        first_line = next(f).strip()
        fields = first_line.split("\t")
        number_of_fields = len(fields)
        
        # Check for the correct number of columns
        if number_of_fields not in {1, 2}:
            raise IndexError(
                f"Input file must have either 1 KO per line or a tab-separated table with 2 columns: "
                f"[id_genome] <tab> [id_ko]. No header expected. Current table contains {number_of_fields} columns."
            )
        
        # Case 1: Single genome (1 column per line)
        if number_of_fields == 1:
            if name is None:
                raise ValueError("If --kos are a single column of KOs then --name must be provided")   
            genome_to_kos[name].add(first_line)
            for line in f:
                line = line.strip()
                if line:
                    genome_to_kos[name].add(line)
                    
        # Case 2: Multiple genomes (2 columns per line)
        elif number_of_fields == 2:
            id_genome, id_ko = fields
            genome_to_kos[id_genome].add(id_ko)
            for line in f:
                line = line.strip()
                if line:
                    id_genome, id_ko = line.split("\t")
                    genome_to_kos[id_genome].add(id_ko)
                    
    # Validate KO identifiers: ensure they are 6 characters long and start with 'K'
    for id_ko in set.union(*genome_to_kos.values()):
        if len(id_ko) != 6 or not id_ko.startswith("K"):
            raise ValueError("Each KO should be 6 characters long and start with 'K'.")
    
    return genome_to_kos

 