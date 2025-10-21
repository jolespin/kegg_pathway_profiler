#!/usr/bin/env python
import sys,os, argparse
from importlib.resources import files as resource_files
from collections import defaultdict
from datetime import datetime
# import pandas as pd 
from tqdm import tqdm

from pyexeggutor import (
    open_file_reader,
    open_file_writer,
    write_pickle, 
    build_logger,
    format_bytes,
)

from kegg_pathway_profiler.pathways import (
#     get_pathway_coverage,
    # pathway_coverage_wrapper,
    Pathway,
)

__program__ = os.path.split(sys.argv[0])[-1]

os.makedirs(resource_files('kegg_pathway_profiler').joinpath('data'), exist_ok=True)
DEFAULT_DATABASE = str(resource_files('kegg_pathway_profiler').joinpath('data/database.pkl.gz'))

# Get the current date and time
now = datetime.now()

# Format the date and time as a string
datetime_string = f"{now.year}.{now.month}.{now.day}"
DATABASE_VERSION = f"KEGG_v{datetime_string}"
    
# http://rest.kegg.jp/list/module
# http://rest.kegg.jp/get/${ID}
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
    parser_io = parser.add_argument_group('Local arguments')
    parser_io.add_argument("-d","--database", type=str,  default=DEFAULT_DATABASE, help = f"path/to/database.pkl[.gz] [Default: {DEFAULT_DATABASE}]")
    parser_io.add_argument("-V", "--database_version", type=str,  default=DATABASE_VERSION, help = f"Database version: Adds version information to the following file: path/to/database.version where .pkl extensions are removed [Default: {DATABASE_VERSION}]")
    parser_io.add_argument("-f", "--force",action="store_true", help = "If file exists, then remove file and update it.")

    parser_local = parser.add_argument_group('Local arguments')
    parser_local.add_argument("-i","--pathway_definitions", type=str, help = "path/to/pathway_definitions.tsv.  [id_pathway]<tab>[definition], No header.")
    parser_local.add_argument("-n","--pathway_names", type=str, help = "path/to/pathway_names.tsv  [id_pathway]<tab>[name], No header.")
    parser_local.add_argument("-c","--pathway_classes", type=str, help = "path/to/pathway_classes.tsv.  [id_pathway]<tab>[class], No header.")
    
    parser_download = parser.add_argument_group('Download arguments')
    parser_download.add_argument("--download", action="store_true",  help = "Download directly from http://rest.kegg.jp/")
    parser_download.add_argument("--intermediate_directory", default="auto", help = "Write the intermediate files from http://rest.kegg.jp/ to a directory.  If 'auto' then download to the directory that contains --database in a subdirectory called `pathway_data`.")
    parser_download.add_argument("--no_intermediate_files", action="store_true",  help = "Don't write intermediate files")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    

    # logger
    logger = build_logger("kegg_pathway_profiler build-pathway-database")

    # Commands
    logger.info(f"Command: {sys.argv}")
    
    # Defaults
    if not opts.download:
        if not opts.pathway_definitions:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
        if not opts.pathway_names:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
        if not opts.pathway_classes:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
    else:
        if opts.pathway_definitions:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
        if opts.pathway_names:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
        if opts.pathway_classes:
            raise ValueError("If --download is not provided then all of the arguments must be provided: -i/--pathway_definitions, -n/--pathway_names, and -c/--pathway_classes")
        
    if os.path.exists(opts.database):
        if not opts.force:
            raise FileExistsError(f"{opts.database} already exists.  To overwrite, please use -f/--force")
        
    # Database version
    database_version_filepath = opts.database
    if database_version_filepath.endswith(".gz"):
        database_version_filepath = database_version_filepath[:-3]
    if database_version_filepath.endswith(".pkl"):
        database_version_filepath = database_version_filepath[:-4]
    if database_version_filepath.endswith(".pickle"):
        database_version_filepath = database_version_filepath[:-7]
    if database_version_filepath.endswith("."):
        database_version_filepath = database_version_filepath[:-1]
    database_version_filepath += ".version"
    database_table_filepath = database_version_filepath[:-8] + ".tsv"
    
    # Intermediate files
    if opts.intermediate_directory == "auto":
        opts.intermediate_directory = os.path.join(os.path.split(opts.database)[0], "pathway_data")
    if not opts.no_intermediate_files:
        os.makedirs(os.path.join(opts.intermediate_directory, "pathways"), exist_ok=True)
    
    # Database
    logger.info(f"Building database version: {opts.database_version}")

    database = defaultdict(dict)
    
    if opts.download:
        from urllib.request import urlopen

        # Fetch the pathway names and read HTML content
        url = "http://rest.kegg.jp/list/module"
        logger.info(f"Fetching KEGG pathway names: {url})")

        with urlopen(url) as response:
            html_content = response.read().decode('utf-8')
            if not opts.no_intermediate_files:
                with open_file_writer(os.path.join(opts.intermediate_directory, "pathway_names.tsv.gz")) as f:
                    print(html_content, file=f)
            for line in html_content.strip().split("\n"):
                line = line.strip()
                if line:
                    id, name = line.split("\t")
                    database[id]["name"] = name
                    
        # Fetch the pathway names and read HTML content
        url = "http://rest.kegg.jp/list/module"
        logger.info(f"Fetching KEGG pathway definitions and classes: {url})")

        # Fetch the HTML content
        for id in tqdm(database, desc=f"Fetching and parsing KEGG"):
            url=f"http://rest.kegg.jp/get/{id}"
            with urlopen(url) as response:
                html_content = response.read().decode('utf-8')
                if not opts.no_intermediate_files:
                    with open_file_writer(os.path.join(opts.intermediate_directory, "pathways", f"{id}.txt.gz")) as f:
                        print(html_content, file=f)
                for line in html_content.strip().split("\n"):
                    line = line.strip()
                    if line:
                        if line.startswith("DEFINITION"):
                            database[id]["definition"] = line[12:]
                        elif line.startswith("CLASS"):
                            database[id]["classes"] = line[12:]
                            
    else:
        # Pathway definitions
        logger.info(f"Reading pathway definitions: {opts.pathway_definitions})")
        
        with open_file_reader(opts.pathway_definitions) as f:
            for line in f:
                line = line.strip()
                if line:
                    id, definition = line.split("\t")
                    database[id]["definition"] = definition
                    
        # Pathway names
        logger.info(f"Reading pathway names: {opts.pathway_names})")
        
        with open_file_reader(opts.pathway_names) as f:
            for line in f:
                line = line.strip()
                if line:
                    id, name = line.split("\t")
                    if id not in database:
                        raise KeyError(f"--pathway_names {opts.pathway_names} contains {id} which is not in --pathway_definitions {opts.pathway_definitions}")
                    database[id]["name"] = name
                    
        # Pathway names
        logger.info(f"Reading pathway classes: {opts.pathway_classes})")
        
        with open_file_reader(opts.pathway_classes) as f:
            for line in f:
                line = line.strip()
                if line:
                    id, classes = line.split("\t")
                    if id not in database:
                        raise KeyError(f"--pathway_classes {opts.pathway_classes} contains {id} which is not in --pathway_definitions {opts.pathway_definitions}")
                    database[id]["classes"] = classes
                    
    # Building
    logger.info(f"Parse pathway definition and building graphs")
    
    for id_pathway, d in tqdm(database.items(), desc=description, unit=" Pathways"):
        # Get attributes
        definition = d["definition"]
        name = d["name"]
        classes = d["classes"]
        # Build and parse pathway
        pathway = Pathway(id=id_pathway, definition=definition, name=name, classes=classes)
        # Store output
        database[id_pathway]["graph"] = pathway.graph_
        database[id_pathway]["ko_to_nodes"] = pathway.ko_to_nodes_
        database[id_pathway]["optional_kos"] = pathway.optional_kos_

   # Write Database
    logger.info(f"Writing database file: {opts.database}")
    write_pickle(database, opts.database)
    
   # Write Database Version
    logger.info(f"Writing database version file: {database_version_filepath}")
    with open_file_writer(database_version_filepath) as f:
        print("VERSION:", opts.database_version, file=f)
        print("CREATED:", now, file=f)
        
   # Write Database KO list
    logger.info(f"Writing database pathway table: {database_table_filepath}")
    with open_file_writer(database_table_filepath) as f:
        for id_pathway, d in tqdm(database.items(), desc=description, unit=" Pathways"):
            for id_ko in d["ko_to_nodes"]:
                print(id_pathway, id_ko, sep="\t", file=f)

    # Summarize database
    size_in_bytes = os.stat(opts.database).st_size
    logger.info(f"Database size: {format_bytes(size_in_bytes)} ({size_in_bytes} bytes)")
    database_kos = set.union(*map(lambda d: set(d["ko_to_nodes"].keys()), database.values()))
    logger.info(f"Number of pathways: {len(database)}")
    logger.info(f"Number of unique KOs: {len(database_kos)}")
    
if __name__ == "__main__":
    main()
    
    

    

