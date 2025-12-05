#!/usr/bin/env python
import sys,os, argparse, warnings, subprocess
from importlib.resources import files as resource_files
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import pandas as pd 
from tqdm import tqdm

from pyexeggutor import (
    read_pickle, 
    write_pickle, 
    build_logger,
    reset_logger,
    format_bytes,
    get_directory_size,
)

from kegg_pathway_profiler.utils import (
    read_kos,
)

from kegg_pathway_profiler.pathways import (
    pathway_coverage_wrapper,
)


__program__ = os.path.split(sys.argv[0])[-1]

DEFAULT_DATABASE = resource_files('kegg_pathway_profiler').joinpath('data/database.pkl.gz')

# Wrapper function for parallel execution
def process_genome(id_genome, genome_kos, database):
    pathway_to_results = pathway_coverage_wrapper(
        evaluation_kos=genome_kos,
        database=database,
        progressbar_description=f"Calculating pathway coverage: {id_genome}",
        progressbar=False,
    )
    
    # Extract coverage data
    coverage_data = {
        id_pathway: results["coverage"] 
        for id_pathway, results in pathway_to_results.items()
    }
    
    return id_genome, pathway_to_results, coverage_data

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
    usage = f"{__program__} "
    epilog = "https://github.com/jolespin/kegg_pathway_profiler"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--kos", type=str, default="stdin", help = "path/to/kos.list[.gz].  Can either be 1 KO per line or a tab-separated table with the following structure: [id_genome]<tab>[id_ko], No header.")
    # parser_io.add_argument("-t","--prevalence_table", type=str,  help = "path/to/prevalence_table.tsv[.gz].  Prevalence table with genomes as rows and KO as columns with counts prevalence values")
    parser_io.add_argument("-n","--name", type=str,  help = "Name of genome. [Default: Filename for --kos]")
    parser_io.add_argument("-o","--output_directory", type=str, default="kegg_pathway_profiler_output", help = "path/to/output_directory/ (e.g., kegg_pathway_profiler_output/]")
    parser_io.add_argument("-d","--database", type=str, default=DEFAULT_DATABASE, help = f"path/to/database.pkl[.gz] [Default: {DEFAULT_DATABASE}]")
    parser_io.add_argument("--index_name", type=str, default="id_genome", help = f"Index name for coverage table (e.g., id_genome, id_genome_cluster, id_contig) [Default: id_genome]")

    # Utilities
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-p","--n_jobs", type=int, default=1,  help = "Number of threads to use.  Use -1 for all available. [Default: 1]")

    # Pathways
    # parser_pathways = parser.add_argument_group('Pathways arguments')
    # parser_pathways.add_argument("-w", "--include_weights", help="Add weights for each KO in output", action='store_true')
    # parser_pathways.add_argument("-p", "--plot_pathways", help="Plot pathway completeness", action='store_true')
    
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # logger
    logger = build_logger("kegg_pathway_profiler profile-pathway-coverage")

    # Commands
    logger.info(f"Command: {sys.argv}")
    

    # Defaults
    if opts.kos == "stdin":
        opts.kos = sys.stdin
        if opts.name is None:
            raise ValueError("If --kos is stdin then --name must be provided")        
        
    # Threads
    if opts.n_jobs == -1:
        logger.info(f"Determining available number of threads")
        from multiprocessing import cpu_count
        opts.n_jobs = cpu_count()
    logger.info(f"Using {opts.n_jobs} threads")

    # Make directory
    logger.info(f"Creating output directory: {opts.output_directory})")
    os.makedirs(opts.output_directory, exist_ok=True)
    
    # Load database
    size_in_bytes = os.stat(opts.database).st_size
    logger.info(f"Database size: {format_bytes(size_in_bytes)} ({size_in_bytes} bytes)")
    logger.info(f"Loading database: {opts.database})")
    database = read_pickle(opts.database)
    database_kos = set.union(*map(lambda d: set(d["ko_to_nodes"].keys()), database.values()))
    logger.info(f"Number of pathways: {len(database)}")
    logger.info(f"Number of unique KOs: {len(database_kos)}")


    # Genome -> KO set
    logger.info(f"Reading query KOs: {opts.kos})")
    genome_to_kos = read_kos(opts.kos, opts.name)

    # Coverage table
    logger.info(f"Calculating pathway coverage")
 
    # Parallel execution
    output_data = {}
    coverage_table = defaultdict(dict)

    with ProcessPoolExecutor(max_workers=opts.n_jobs) as executor:
        # Submit all jobs
        futures = {
            executor.submit(process_genome, id_genome, genome_kos, database): id_genome
            for id_genome, genome_kos in genome_to_kos.items()
        }
        
        # Collect results with progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing genomes"):
            id_genome, pathway_to_results, coverage_data = future.result()
            output_data[id_genome] = pathway_to_results
            coverage_table[id_genome] = coverage_data

    # Coverage table
    df_coverage_table = pd.DataFrame(coverage_table).T.fillna(0.0)
    df_coverage_table = df_coverage_table.loc[:,sorted(df_coverage_table.columns, key=lambda id_pathway:int(id_pathway[1:]))]
    df_coverage_table.index.name = opts.index_name

    output_filepath = os.path.join(opts.output_directory, "pathway_coverage.tsv.gz")
    logger.info(f"Writing pathway coverage table: {output_filepath}")
    df_coverage_table.to_csv(output_filepath, sep="\t")
    
    # Pathway outputs
    output_filepath = os.path.join(opts.output_directory, "pathway_output.pkl.gz")
    logger.info(f"Writing pathway output pickle: {output_filepath}")
    write_pickle(dict(output_data), output_filepath)
    
if __name__ == "__main__":
    main()
    
    

    

