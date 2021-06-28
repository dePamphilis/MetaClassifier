#!/usr/bin/env python3

"""
metaclassifier.process_reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions for merging paired-end (PE) and converting FASTQ to FASTA

A collection of functions for merging overlapping PE reads when the DNA fragment is shorter
than two times the read length, and converting the resulting merged FASTQ to FASTA format.
"""

import os
import sys
import time
import subprocess
from datetime import date
from . import meta_common as utils 

def merge_pairs(sample_input, output_dir, merger, threads):
    """Merges overlapping paired-end (PE) sample FASTQ reads

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    output_dir : str
        output directory
    merger : str
        The PEAR, the PE read merger executable
    threads : str 
        The number of threads for PEAR, the PE read merger

    Returns
    -------
        None
    """

    print(f"{utils.get_localtime()} - - Merging paired-end (PE) sample read dataset(s)...\n")

    # get paired-end sample read file and merge with PEAR
    merge_dir = f"{output_dir}/merge_dir"
    os.mkdir(merge_dir,  mode=0o775)
    with open(f"{output_dir}/sample.tsv", 'w') as output_file:
        for sample in utils.parse_config_file(sample_input):
            merged_fastq = f"{merge_dir}/{sample[0]}"
            log_file = f"{merge_dir}/{sample[0]}_pear.log"
            with open(log_file, 'w') as log:
                # TODO: replace this with something better
                subprocess.run([merger, '-f', sample[1], '-r', sample[2], '-o', merged_fastq, '-j', threads], stdout=log, check=True)
            output_file.write(f"{sample[0]}\t{merge_dir}/{sample[0]}.assembled.fastq\n")

def get_fasta(sample_input, output_dir, converter):
    """Converts merged paired-end (PE) reads from FASTQ to FASTA 

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    output_dir : str
        output directory
    converter : str
        The seqtk, the sequence processing tool executable

    Returns
    -------
        None
    """
    print(f"{utils.get_localtime()} - - Converting merged paired-end (PE) sample read dataset(s) from FASTQ to FASTA...\n")

    # convert merged paired-end fastq sample reads files to fasta format
    fasta_dir = f"{output_dir}/fasta_dir"
    os.mkdir(fasta_dir,  mode=0o775)
    for sample in utils.parse_config_file(sample_input, 2):
        fasta_file = f"{fasta_dir}/{sample[0]}.fasta"
        with open(fasta_file, 'w') as fasta:
            subprocess.run([converter, 'seq', '-A', sample[1]], stdout=fasta, check=True)

def main():
    t0 = time.time()
    args = utils.read_parameters("process")

    print(f"{utils.get_localtime()} - Starting read processing...\n")
    if args.frag_type == "single" and args.merge == True:
        raise Exception(f"single-end read fragments cannot be merged!")
    if args.frag_type == "paired" and args.merge == False:
        raise Exception(f"Overlapping paired-end read fragments need to be merged!") 

    if not args.output_dir:
        OUTPUT_DIR = "reads_output_"+str(date.today())
        os.mkdir(OUTPUT_DIR,  mode=0o775)
    else:
        OUTPUT_DIR = args.output_dir
        if not os.path.isdir(OUTPUT_DIR):
            os.mkdir(OUTPUT_DIR,  mode=0o775)

    if args.merge == True:
        merge_pairs(args.SAMPLE_FILE, OUTPUT_DIR, args.pear_merger, args.threads)
        sample_input = f"{OUTPUT_DIR}/sample.tsv"
        get_fasta(sample_input, OUTPUT_DIR, args.seqtk_converter)
    else:
        get_fasta(args.SAMPLE_FILE, OUTPUT_DIR, args.seqtk_converter) 

    t1 = time.time()
    print(f"{utils.get_localtime()} - Completed read processing...")
    print(f"Total elapsed time {int(t1 - t0)}\n")
    sys.exit(0)

if __name__ == '__main__':
    main()
