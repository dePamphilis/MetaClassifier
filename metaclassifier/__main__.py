#!/usr/bin/env python3

"""
metaclassifier.__main___
~~~~~~~~~~~~~~~~~~~~~~~~

Function quantify the taxonomic floral composition of honey samples read datasets

A wrapper module the "process_reads" and "classify_reads" modules "for merging paired-end (PE)
and converting FASTQ to FASTA" and "for classifying and quantifying sample reads into taxonomy
groups" respectively.
"""


import os
import sys
import time
from . import meta_common as utils
from . import process_reads, classify_reads


# improvements
# Allow multiple file types for input/output (csv,tsv,json,etc...)
# Validate input file before running program, check file is valid, and check all files exist

def run_metaclassifier(sample_input, db_dir, config_file, output_dir, frag_type, merge, merger,  
                    converter, aligner, tax_class, min_proportion, max_markers, threads):
    """Merges overlapping paired-end (PE) sample FASTQ reads

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    db_dir : str 
        Input marker database directory
    config_file : str
        Input tab-delimited file specifying marker name, and its
        corresponding VSEARCH search parameters
    output_dir : str
        output directory
    frag_type : str
        The sequence fragment type in the input sample file
        (either "paired" or "single")
    merge : boolean
        Either True or False
    merger : str
        The PEAR, the PE read merger executable
    converter : str
        The seqtk, the sequence processing tool executable
    aligner : str
        VSEARCH, the sequence analysis tool executable
    tax_class : str
        Taxonomy class for quantify taxon level marker read abundance
    min_proportion : float
        Minimum taxon read proportion allowed to retain a sample taxon
    max_markers = int
        Maximum missing markers allowed to retain a sample taxon    
    threads : str 
        The number of threads for PEAR, the PE read merger

    Returns
    -------
        None
    """

    print(f"\n===========================================================")
    print(f"MetaClassifier version {utils.__version__} (release date: {utils.__date__})")
    print(f"===========================================================\n")
    print(f"{utils.get_localtime()} - Starting MetaClassifier...\n")

    if frag_type == "single" and merge == True:
        raise Exception(f"single-end read frangments cannot be merged!")

    if frag_type == "paired" and merge == False:
        raise Exception(f"Overlapping paired-end read frangments need to be merged!")

    if min_proportion >= 0.01:
        raise Exception(f"The minimum taxon read proportion allowed to retain a sample taxon cannot exceed 0.1%!")

    if merge == True:
        process_reads.merge_pairs(sample_input, output_dir, merger, threads)
        sample_input = f"{output_dir}/sample.tsv"
        process_reads.get_fasta(sample_input, output_dir, converter)
    else:
        process_reads.get_fasta(sample_input, output_dir, converter)

    fasta_dir = f"{output_dir}/fasta_dir"
    classify_reads.search_markers(fasta_dir, db_dir, config_file, aligner, tax_class, min_proportion, max_markers, threads)

    print(f"{utils.get_localtime()} - Completed MetaClassifier...")

def main():
    t0 = time.time()
    
    args = utils.read_parameters("main")
    if not args.output_dir:
        OUTPUT_DIR = os.path.splitext(os.path.basename(args.SAMPLE_FILE))[0]
        os.mkdir(OUTPUT_DIR,  mode=0o775)
    else:
        OUTPUT_DIR = args.output_dir
    
    run_metaclassifier(args.SAMPLE_FILE, args.DB_DIR, args.CONFIG_FILE, OUTPUT_DIR, args.frag_type, args.merge, args.pear_merger, 
                    args.seqtk_converter, args.vsearch_aligner, args.tax_class, args.min_proportion, args.max_markers, args.threads)

    t1 = time.time()
    print(f"Total elapsed time {int(t1 - t0)}\n")
    sys.exit(0)

if __name__ == '__main__':
    main()
    
