#!/usr/bin/env python
import sys,os, argparse
from collections import defaultdict
from tqdm import tqdm
from pyexeggutor import (
    open_file_reader,
    open_file_writer,
)

__program__ = os.path.split(sys.argv[0])[-1]


def main(args=None):
    # Options
    # =======
    # Path info
    python_executable = sys.executable
    bin_directory = "/".join(python_executable.split("/")[:-1])
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, sys.version.split(" ")[0], python_executable, script_filename)
    usage = f"{__program__} -i <pykofamsearch_output.tsv> -o <kos.list>"
    epilog = "https://github.com/jolespin/kegg_pathway_profiler"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i", "--input", type=str, default="stdin", help = "path/to/pykofamsearch_output.tsv [Default: stdin]")
    parser.add_argument("-o", "--output", type=str, default="stdout", help = "Either 1) list of KOs; or 2) [id_genome]<tab>[id_ko] if --name/--identifier_mapping is provided. No header. [Default: stdout]")
    parser.add_argument("-m", "--identifier_mapping", type=str, help = "Identifier mapping [id_protein]<tab>[id_genome] or [id_protein]<tab>[id_contig]<tab>[id_genome], No header. Cannot be used with --name.")
    parser.add_argument("-f", "--identifier_mapping_format", type=int,  choices={1,2,3}, help = "1: --identifier_mapping <id_genome>, 2: path/to/tsv [id_protein]<tab>[id_genome], or 3: path/to/tsv [id_protein]<tab>[id_contig]<tab>[id_genome]")
    parser.add_argument("-r", "--reformatted", action="store_true", help = "PyKOfamSearch output is reformatted")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # I/O
    if opts.input == "stdin":
        f_in = sys.stdin
    else:
        f_in = open_file_reader(opts.input)
        next(f_in)
    if opts.output == "stdout":
        f_out = sys.stdout
    else:
        f_out = open_file_writer(opts.output)
    
    if opts.identifier_mapping:
        if not opts.identifier_mapping_format:
            raise ValueError("If -m/--identifier_mapping is provided then you must also provide -f/--identifier_mapping_format")
        
        # Protein --> Genome mapping
        protein_to_genome = dict()
        if opts.identifier_mapping_format > 1:
            # Read identifier mapping
            # [id_protein]<tab>[id_genome]
            if opts.identifier_mapping_format == 2:
                with open_file_reader(opts.identifier_mapping) as f:
                    for line in tqdm(f):
                        line = line.strip()
                        if line:
                            id_protein, id_genome = line.split("\t")
                            protein_to_genome[id_protein] = id_genome
                            
            # [id_protein]<tab>[id_contig]<tab>[id_genome]
            elif opts.identifier_mapping_format == 3:
                with open_file_reader(opts.identifier_mapping) as f:
                    for line in tqdm(f):
                        line = line.strip()
                        if line:
                            id_protein, id_contig, id_genome = line.split("\t")
                            protein_to_genome[id_protein] = id_genome 
                            

        # Read PyKOfamSearch
        genome_to_kos = defaultdict(set)
        if opts.identifier_mapping_format > 1:
            if opts.reformatted:
                for line in tqdm(f_in):
                    line = line.strip()
                    if line:
                        id_protein, number_of_hits, ids, *_ = line.split("\t")
                        ids = eval(ids)
                        id_genome = protein_to_genome[id_protein]
                        genome_to_kos[id_genome] |= set(ids)
            else:
                for line in tqdm(f_in):
                    line = line.strip()
                    if line:
                        id_protein, id_ko, *_ = line.split("\t")
                        id_genome = protein_to_genome[id_protein]
                        genome_to_kos[id_genome].add(id_ko)
        else:
            if opts.reformatted:
                for line in tqdm(f_in):
                    line = line.strip()
                    if line:
                        id_protein, number_of_hits, ids, *_ = line.split("\t")
                        ids = eval(ids)
                        id_genome = opts.identifier_mapping
                        genome_to_kos[id_genome] |= set(ids)
            else:
                for line in tqdm(f_in):
                    line = line.strip()
                    if line:
                        id_protein, id_ko, *_ = line.split("\t")
                        id_genome = opts.identifier_mapping
                        genome_to_kos[id_genome].add(id_ko)
                        
        # Output
        for id_genome, kos in tqdm(genome_to_kos.items()):
            for id_ko in kos:
                print(id_genome, id_ko, sep="\t", file=f_out)
                
    else:
        # No identifier mapping - just extract all unique KO IDs
        all_ids = set()
        if opts.reformatted:
            for line in tqdm(f_in):
                line = line.strip()
                if line:
                    id_protein, number_of_hits, ids, *_ = line.split("\t")
                    ids = eval(ids)
                    all_ids |= set(ids)
        else:
            for line in tqdm(f_in):
                line = line.strip()
                if line:
                    id_protein, id_ko, *_ = line.split("\t")
                    all_ids.add(id_ko)
        all_ids = sorted(all_ids)
    
        # Output
        for id_ko in tqdm(all_ids):
            print(id_ko, file=f_out)
    
    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()
    
    
if __name__ == "__main__":
    main()