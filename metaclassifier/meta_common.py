#!/usr/bin/env python3
# meta_common.py

"""
metaclassifier.meta_common
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Common Functions used across the metaclassifier project

A collection of functions used to help the other metaclassifier functions run, and limit the amount of
rewritten code.
"""

import argparse
import time
import csv
import sys
from urllib import request
import tarfile
from os import path, listdir, remove, mkdir, rename, rmdir
from datetime import date
from . import __version__ as version

__author__ = ('Eric Wafula (ewafula@gmail.com)')
__version__ = version
__date__ = '05 March 2021'

OUTPUT_DIR = ""
FRAGMENT_CHOICES = ['paired', 'single']
TAX_CLASS_CHOICES = ["order", "family", "genus", "species"]
MIN_PROPORTION = 0.00
MAX_MARKERS = 0
PEAR = "pear"
SEQTK = "seqtk"
VSEARCH = "vsearch"
THREADS = "2"

def read_parameters(callingClass:str="main"):
    """Read CLI parameters

    Parameters
    ----------
        None

    Returns
    -------
    args : 
        SAMPLE_FILE: str
        FASTA_DIR: str
        DB_DIR: str
        CONFIG_FILE: str
        output_dir: str
        frag_type: str
        merge: bool
        tax_class: str
        min_proportion:float
        max_markers: float
        pear_merger: str
        seqtk_converter: str
        vsearch_aligner: str
        threads: str
    """

    if(callingClass not in ["main","classify","process"]):
        print(callingClass + " is not an acceptable value, please use:\n[\"main\",\"classify\",\"process\"]")

    p = argparse.ArgumentParser(description=("The metaclassifier.py uses DNA metabarcoding sequence reads and"        
                                             " a database of marker sequences, including\ncorresponding taxonomy"
                                             " classes to identify and quantify the floral composition of honey"),
                                formatter_class=argparse.RawTextHelpFormatter)
    
    if(callingClass in ["main","process"]):
        p.add_argument('SAMPLE_FILE', type=str, default=None,
                        help="Input tab-delimited file specifying sample names, file names for forward paired-end\n"
                            "reads, and file names for reverse paired-end reads (full file path if not in working directory)\n"
                            "The second file not required for single-end frangments\n\n")
    
    if(callingClass == "classify"):
        p.add_argument('FASTA_DIR', type=str, default=None,
                        help="Input sample read data fasta file directory - either merged paired-end PE or\n"
                            "single-end fasta files\n\n")

    if(callingClass in ["main","classify"]):
        p.add_argument('DB_DIR', type=str, default=None,
                        help="Input marker database directory with sequence fasta and corresponding taxonomy lineage\n"
                            "files for each marker\n\n")
        p.add_argument('CONFIG_FILE', type=str, default=None,
                        help="Input tab-delimited file specifying marker name, and its corresponding VSEARCH's\n"
                            "usearch_global function minimum query coverage (i.e. 0.8 for 80%%) and minimun sequence\n"
                            "identity (i.e. 0.95 for 95%%) for each search marker (provide the full file path of the\n"
                            "VSEARCH settings configurations if is not in the working directory)\n\n")
    if(callingClass in ["main","process"]):
        p.add_argument('-o', '--output_dir', type=str, default=OUTPUT_DIR,
                        help="Specify output directory name, otherwise it will automatically be created using the\n"
                            "input sample table file name\n\n")   
        p.add_argument('-f', '--frag_type', type=str, default='paired', choices=FRAGMENT_CHOICES,
                        help="Specify the sequence fragment type in the input sample file, available options are:\n"
                            "paired: single-end read fragments (default)\n"
                            "single: paired-end read fragments\n\n")
        p.add_argument('-m', '--merge', action='store_true', 
                        help="Merge overlapping paired-end reads (default: False)\n\n")
    if(callingClass in ["main","classify"]):
        p.add_argument('-c', '--tax_class', type=str, default='genus', choices=TAX_CLASS_CHOICES,
                        help="Taxonomy class for quantify taxon level marker read abundance (default: genus)\n\n")
        p.add_argument('-p', '--min_proportion', type=float, default=MIN_PROPORTION,
                        help="Minimum taxon read proportion allowed to retain a sample taxon, allowed proportion,\n"
                            "ranges from 0.00 to 0.01 (default = 0.00)\n\n")
        p.add_argument('-i', '--max_markers', type=int, default=MAX_MARKERS,
                        help="Maximum missing markers allowed to retain a sample taxon (default = 0)\n\n")
    if(callingClass in ["main","process"]):
        p.add_argument('-r', '--pear_merger', type=str, default=PEAR,
                        help="Path to PEAR, the paired-end read merger if not in environmental variables (ENV)\n"
                            "(default: read from ENV)\n\n")
        p.add_argument('-s', '--seqtk_converter', type=str, default=SEQTK,
                        help="Path to seqtk, the sequence processing tool if not in environmental variables (ENV)\n"
                            "(default: read from ENV)\n\n")
    if(callingClass in ["main", "classify"]):
        p.add_argument('-a', '--vsearch_aligner', type=str, default=VSEARCH,
                        help="Path to VSEARCH, the sequence analysis tool if not in environmental variables (ENV)\n"
                            "(default: read from ENV)\n\n")                 
    p.add_argument('-t', '--threads', type=str, default=THREADS,
                    help="Specify the number of threads to use (default: 2)\n\n")
    p.add_argument('-v', '--version', action='version',
                   version=f"metaclassifier.py version {__version__} ({__date__})",
                   help="Print the current metaclassifier.py version and exit\n\n")
    if(not len(sys.argv) > 1):
        p.print_help()
        sys.exit(0)

    output_args = p.parse_args()
    if(callingClass in ["main","process"]):
        if path.exists(output_args.DB_DIR):
            db_files = listdir(output_args.DB_DIR)
            if any(".fa" in file for file in db_files):
                return output_args
        print("No db files found, downloading defaults")
        output_args.DB_DIR = output_args.DB_DIR+"_defaults"
        download_db(output_args.DB_DIR)


    return output_args

def parse_config_file(configFileName:str, rowLen:int = 3, areFiles:bool = True):
    """Preparse input file check for errors, and return object
       if areFiles is true, assums 1 -> rowLen-1 are files and checks they exist
    
    Parameters:
        configFileName: str
            Path to input file
        rowLen: int default 3
            number of expected rows in file
        areFiles: bool default True
            Bool to decide if we need to check for file validity

    Output:
        fileContents: list of lists
    """
    fileContents = []

    with open(configFileName, 'r') as inputFile:
        delimiter = '\t'
        if (configFileName.split('.'))[-1] == 'csv':
            delimiter = ','
        rows = csv.reader(inputFile,delimiter=delimiter)
        i = 0
        for row in rows:
            i+=1
            if len(row) != rowLen:
                raise ValueError("In Row" + str(i) + "the number of elements does no equal 3")
            row = [x.strip() for x in row]
            if areFiles:
                for x in range(1,rowLen):
                    if not path.exists(row[x]):
                        raise ValueError("File at: " + row[x] + " does not exist")
            fileContents.append(row)

    return fileContents

def download_db(db_path:str, url:str="http://bigdata.bx.psu.edu/MetaClassifier_databases/MetabarcodeDBsV2.tar.gz"):
    """

    Parameters
    ----------
        db_path: str
            path to db files
        url: str
            url of db files

    Returns
    -------
        None
    """

    if not path.exists(db_path):
        mkdir(db_path)
    request.urlretrieve(url,"MetabarcodeDBsV2.tar.gz")
    tar = tarfile.open("./MetabarcodeDBsV2.tar.gz", "r:gz")
    tar.extractall(db_path)
    dir_files = listdir(db_path)
    if len(dir_files) == 1 and path.isdir(db_path+"/"+dir_files[0]):
        full_path = db_path+"/"+dir_files[0]
        for db_file in listdir(full_path):
            rename(full_path+"/"+db_file,db_path+"/"+db_file)
        rmdir(full_path)

    remove("./MetabarcodeDBsV2.tar.gz")


def get_localtime():
    """Get Current Local Time

    Parameters:
        None

    Output:
        current_time: str

    """
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return current_time