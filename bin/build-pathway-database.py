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

EBI_REPO = "EBI-Metagenomics/kegg-pathways-completeness-tool"
EBI_RAW_BASE = f"https://raw.githubusercontent.com/{EBI_REPO}"
EBI_API_BASE = f"https://api.github.com/repos/{EBI_REPO}"
EBI_TARBALL_BASE = f"https://github.com/{EBI_REPO}/archive/refs/tags"
EBI_DATA_SUBPATH = "kegg_pathways_completeness/pathways_data"

def _parse_ebi_modules_table(content, logger):
    database = defaultdict(dict)
    lines = content.strip().split('\n')
    logger.info(f"Parsing unified modules table ({len(lines) - 1} entries)")
    for line in lines[1:]:
        line = line.strip()
        if line:
            fields = line.split('\t')
            if len(fields) != 4:
                raise ValueError(f"Expected 4 tab-separated fields, got {len(fields)}: {line[:80]}...")
            module_id, definition, name, classes = fields
            database[module_id]["definition"] = definition
            database[module_id]["name"] = name
            database[module_id]["classes"] = classes
    return database

def _parse_ebi_legacy_file(content):
    result = {}
    for line in content.strip().split('\n'):
        line = line.strip()
        if line:
            module_id, value = line.split(':', 1)
            result[module_id] = value
    return result

def _fetch_ebi_legacy_from_branch(branch_name, logger):
    from urllib.request import urlopen
    from urllib.error import HTTPError

    database = defaultdict(dict)
    file_map = {
        "all_pathways.txt": "definition",
        "all_pathways_names.txt": "name",
        "all_pathways_class.txt": "classes",
    }
    for filename, field in file_map.items():
        url = f"{EBI_RAW_BASE}/{branch_name}/{EBI_DATA_SUBPATH}/{filename}"
        logger.info(f"Fetching EBI legacy file: {url}")
        with urlopen(url) as response:
            content = response.read().decode('utf-8')
        parsed = _parse_ebi_legacy_file(content)
        for module_id, value in parsed.items():
            database[module_id][field] = value
    return database

def _parse_ebi_legacy_from_tarball(tar, data_prefix, logger):
    database = defaultdict(dict)
    file_map = {
        "all_pathways.txt": "definition",
        "all_pathways_names.txt": "name",
        "all_pathways_class.txt": "classes",
    }
    for filename, field in file_map.items():
        member_path = f"{data_prefix}/{filename}"
        logger.info(f"Extracting EBI legacy file: {member_path}")
        f = tar.extractfile(member_path)
        content = f.read().decode('utf-8')
        parsed = _parse_ebi_legacy_file(content)
        for module_id, value in parsed.items():
            database[module_id][field] = value
    return database

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
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-d","--database", type=str,  default=DEFAULT_DATABASE, help = f"path/to/database.pkl[.gz] [Default: {DEFAULT_DATABASE}]")
    parser_io.add_argument("-V", "--database_version", type=str,  default=DATABASE_VERSION, help = f"Database version: Adds version information to the following file: path/to/database.version where .pkl extensions are removed [Default: {DATABASE_VERSION}]")
    parser_io.add_argument("-f", "--force",action="store_true", help = "If file exists, then remove file and update it.")
    parser_io.add_argument("--intermediate_directory", default="auto", help = "Write intermediate files to a directory.  If 'auto' then use the directory that contains --database in a subdirectory called `pathway_data`.")
    parser_io.add_argument("--no_intermediate_files", action="store_true",  help = "Don't write intermediate files")

    parser_ebi = parser.add_argument_group('EBI repository arguments (recommended)')
    parser_ebi.add_argument("--ebi", type=str, default=None,
        help = "Download from EBI kegg-pathways-completeness-tool GitHub repository (recommended).\n"
               "Accepts: 'latest' (latest release), a release tag (e.g., '1.4.3'),\n"
               "or 'branch:<name>' (e.g., 'branch:master').\n"
               "Examples: --ebi latest, --ebi 1.4.3, --ebi branch:master")

    parser_local = parser.add_argument_group('Local arguments')
    parser_local.add_argument("-i","--pathway_definitions", type=str, help = "path/to/pathway_definitions.tsv.  [id_pathway]<tab>[definition], No header.")
    parser_local.add_argument("-n","--pathway_names", type=str, help = "path/to/pathway_names.tsv  [id_pathway]<tab>[name], No header.")
    parser_local.add_argument("-c","--pathway_classes", type=str, help = "path/to/pathway_classes.tsv.  [id_pathway]<tab>[class], No header.")

    parser_download = parser.add_argument_group('Download arguments')
    parser_download.add_argument("--download", action="store_true",  help = "Download directly from http://rest.kegg.jp/")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename


    # logger
    logger = build_logger("kegg_pathway_profiler build-pathway-database")

    # Commands
    logger.info(f"Command: {sys.argv}")

    # Validate mutual exclusivity of modes
    has_local = any([opts.pathway_definitions, opts.pathway_names, opts.pathway_classes])
    modes = []
    if opts.download:
        modes.append("--download")
    if opts.ebi is not None:
        modes.append("--ebi")
    if has_local:
        modes.append("local files (-i/-n/-c)")

    if len(modes) > 1:
        raise ValueError(f"Mutually exclusive modes specified: {', '.join(modes)}. "
                         "Use exactly one of: --download, --ebi, or -i/-n/-c local files.")

    if len(modes) == 0:
        raise ValueError("Must specify one of: --download, --ebi, or all three local file "
                         "arguments (-i/--pathway_definitions, -n/--pathway_names, -c/--pathway_classes)")

    if has_local:
        if not opts.pathway_definitions:
            raise ValueError("Local file mode requires all three arguments: -i/--pathway_definitions, -n/--pathway_names, -c/--pathway_classes")
        if not opts.pathway_names:
            raise ValueError("Local file mode requires all three arguments: -i/--pathway_definitions, -n/--pathway_names, -c/--pathway_classes")
        if not opts.pathway_classes:
            raise ValueError("Local file mode requires all three arguments: -i/--pathway_definitions, -n/--pathway_names, -c/--pathway_classes")

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

    elif opts.ebi is not None:
        from urllib.request import urlopen
        from urllib.error import HTTPError
        import tarfile
        import io
        import json

        ebi_spec = opts.ebi

        if ebi_spec.startswith("branch:"):
            branch_name = ebi_spec[len("branch:"):]
            fetch_mode = "branch"
        elif ebi_spec == "latest":
            url = f"{EBI_API_BASE}/releases/latest"
            logger.info(f"Fetching latest EBI release tag: {url}")
            try:
                with urlopen(url) as response:
                    release_info = json.loads(response.read().decode('utf-8'))
                    tag = release_info["tag_name"]
            except HTTPError as e:
                raise ConnectionError(f"Failed to fetch latest EBI release info from {url}: {e}") from e
            logger.info(f"Latest EBI release: {tag}")
            fetch_mode = "release"
        else:
            tag = ebi_spec
            fetch_mode = "release"

        if fetch_mode == "branch":
            modules_table_url = f"{EBI_RAW_BASE}/{branch_name}/{EBI_DATA_SUBPATH}/modules_table.tsv"
            logger.info(f"Attempting to fetch unified modules table: {modules_table_url}")
            try:
                with urlopen(modules_table_url) as response:
                    content = response.read().decode('utf-8')
                database = _parse_ebi_modules_table(content, logger)
            except HTTPError:
                logger.info("Unified table not found, trying legacy 3-file format")
                database = _fetch_ebi_legacy_from_branch(branch_name, logger)

        elif fetch_mode == "release":
            tarball_url = f"{EBI_TARBALL_BASE}/{tag}.tar.gz"
            logger.info(f"Downloading EBI release tarball: {tarball_url}")
            try:
                with urlopen(tarball_url) as response:
                    tarball_bytes = response.read()
            except HTTPError as e:
                if e.code == 404:
                    raise ValueError(
                        f"EBI release '{tag}' not found. "
                        f"Check available releases at https://github.com/{EBI_REPO}/releases"
                    ) from e
                raise

            with tarfile.open(fileobj=io.BytesIO(tarball_bytes), mode='r:gz') as tar:
                members = tar.getmembers()
                top_dir = members[0].name.split('/')[0]
                data_prefix = f"{top_dir}/{EBI_DATA_SUBPATH}"
                member_names = {m.name for m in members}

                if f"{data_prefix}/modules_table.tsv" in member_names:
                    f = tar.extractfile(f"{data_prefix}/modules_table.tsv")
                    content = f.read().decode('utf-8')
                    database = _parse_ebi_modules_table(content, logger)
                elif f"{data_prefix}/all_pathways.txt" in member_names:
                    database = _parse_ebi_legacy_from_tarball(tar, data_prefix, logger)
                else:
                    raise FileNotFoundError(
                        f"EBI release {tag} does not contain expected data files in {EBI_DATA_SUBPATH}/")

        # Save intermediate files
        if not opts.no_intermediate_files:
            logger.info(f"Writing intermediate files to: {opts.intermediate_directory}")
            with open_file_writer(os.path.join(opts.intermediate_directory, "pathway_definitions.tsv.gz")) as f:
                for id_pathway in sorted(database):
                    print(f"{id_pathway}\t{database[id_pathway]['definition']}", file=f)
            with open_file_writer(os.path.join(opts.intermediate_directory, "pathway_names.tsv.gz")) as f:
                for id_pathway in sorted(database):
                    print(f"{id_pathway}\t{database[id_pathway]['name']}", file=f)
            with open_file_writer(os.path.join(opts.intermediate_directory, "pathway_classes.tsv.gz")) as f:
                for id_pathway in sorted(database):
                    print(f"{id_pathway}\t{database[id_pathway]['classes']}", file=f)

        # Auto-set database version if user didn't override
        if opts.database_version == DATABASE_VERSION:
            if fetch_mode == "branch":
                opts.database_version = f"ebi-{branch_name}_v{datetime_string}"
            else:
                opts.database_version = f"ebi_v{tag}"

        logger.info(f"Loaded {len(database)} pathways from EBI repository")

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

    # Validate database entries
    ids_to_remove = []
    for id_pathway, d in database.items():
        missing_fields = [f for f in ["definition", "name", "classes"] if f not in d]
        if missing_fields:
            logger.warning(f"Module {id_pathway} is missing fields: {', '.join(missing_fields)}. Skipping.")
            ids_to_remove.append(id_pathway)
    for id_pathway in ids_to_remove:
        del database[id_pathway]
    if not database:
        raise ValueError("No valid pathway entries found")

    # Building
    logger.info(f"Building database version: {opts.database_version}")
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
