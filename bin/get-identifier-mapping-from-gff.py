#!/usr/bin/env python
import sys,os, argparse
from collections import defaultdict
from tqdm import tqdm
from pyexeggutor import (
    open_file_reader,
    open_file_writer,
    parse_attribute_from_gff,
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
    usage = f"{__program__} "
    epilog = "https://github.com/jolespin/kegg_pathway_profiler"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i", "--input", type=str, default="stdin", help = "path/to/file.gff[.gz] [Default: stdin]")
    parser.add_argument("-o", "--output", type=str, default="stdout", help = "Either 1) [id_protein]<tab>[id_contig] or 2) [id_protein]<tab>[id_contig]<tab>[id_genome] if --name is provided. [Default: stdout]")
    parser.add_argument("-n", "--genome", type=str,  help = "ID of genome")
    parser.add_argument("-a", "--attribute_key", type=str,  help = "Attribute key to parse from GFF file (e.g., ID, Name, gene_id, etc.) [Default: ID]", default="ID")
    parser.add_argument("-t", "--feature_type", type=str,  help = "Feature type to parse from GFF file (e.g., gene, CDS, mRNA, etc.) [Default: CDS]", default="CDS")
    parser.add_argument("--header", action="store_true", help = "Output has header")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # I/O
    if opts.input == "stdin":
        f_in = sys.stdin
    else:
        f_in = open_file_reader(opts.input)

    if opts.output == "stdout":
        f_out = sys.stdout
    else:
        f_out = open_file_writer(opts.output)


    if opts.genome:
        if opts.header:
            print("id_protein", "id_contig", "id_genome", sep="\t", file=f_out)
        for id_contig, id_protein in parse_attribute_from_gff(f_in, attribute_key=opts.attribute_key, feature_type=opts.feature_type):
            print(id_protein, id_contig, opts.genome, sep="\t", file=f_out)
    else:
        if opts.header:
            print("id_protein", "id_contig", sep="\t", file=f_out)
        for id_contig, id_protein in parse_attribute_from_gff(f_in, attribute_key=opts.attribute_key, feature_type=opts.feature_type):
            print(id_protein, id_contig, sep="\t", file=f_out)
            
    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()
    
    
if __name__ == "__main__":
    main()