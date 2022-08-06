#!/usr/bin/env python3

"""
metaclassifier.classify_reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions for classifying and quantifying sample reads into taxonomy groups

A collection of functions to searching for marker sequences in mult-fasta short reads sample
datasets, identifying sample taxonomic composition, and quantifying taxon read abundance.
"""

import os
import sys
import time
import subprocess
import argparse
import pandas as pd
from functools import reduce
from . import meta_common as utils 

def search_markers(fasta_dir, db_dir, config_file, aligner, tax_class, min_proportion, max_markers, threads):
     """Searches for marker sequences in mult-fasta short reads sample datasets

     Parameters
     ----------
     fasta_dir : str
          Input sample read data fasta file directory
     db_dir : str 
          Input marker database directory
     config_file : str
          Input tab-delimited file specifying marker name, and its
          corresponding VSEARCH search parameters
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

     print(f"{utils.get_localtime()} - - Searching for markers in the sample read dataset(s)...\n")

     # get sample name and create sample out directory
     for file in os.listdir(fasta_dir):
          sample_fasta = f"{fasta_dir}/{file}"
          sample = os.path.splitext(file)[0]
          print(f"{utils.get_localtime()} - - - Searching sample {sample}...\n")
          fasta_dir_path_components = os.path.abspath(fasta_dir).split(os.sep)
          output_dir = "/".join(fasta_dir_path_components[:-1])
          sample_dir = f"{output_dir}/{sample}"
          os.mkdir(sample_dir,  mode=0o775)
          

          # search for marker with vsearch aligner
          marker_dfs = []  
          for params in utils.parse_config_file(config_file,3,False):
               marker = params[0]
               seq_identity = params[1]
               query_coverage = params[2]
               marker_fasta = f"{db_dir}/{params[0]}.fa"
               marker_search = f"{sample_dir}/{params[0]}_blast6out.tsv"
               log_file = f"{sample_dir}/{params[0]}_vsearch.log"
               subprocess.run([aligner, '--usearch_global', sample_fasta, '--db', marker_fasta, '--blast6out', marker_search, '--id', seq_identity, '--query_cov', query_coverage, '--top_hits_only', '--log', log_file, '--quiet', '--threads', threads], capture_output=True, check=True)
               
               # get taxon read proportions
               marker_lineage = f"{db_dir}/{marker}.tax"
               marker_dfs.append(get_read_proportion(marker_search, marker_lineage, tax_class, min_proportion, marker))
                    
          # get rescaled median taxon read proportions and write final sample results (markers and taxa) to file
          sample_df = get_rescaled_median_proportion(marker_dfs, tax_class, max_markers)
          sample_results = f"{sample_dir}/{sample}_rescaled_propotions.tsv"
          sample_df.to_csv(sample_results, sep="\t", index=False, float_format="%.6f", encoding="utf-8")

def get_read_proportion(align_resluts, marker_lineage, tax_class, min_proportion, marker):
     """Identifies sample taxonomic composition and quantifies class read abundance

     Parameters
     ----------
     align_resluts : str
          VSEARCH blast6out sequence search results
     marker_lineage : str 
          the marker lineage table from the input marker database directory
     tax_class : str
          Taxonomy class for quantify taxon level marker read abundance
     min_proportion : float
          Minimum taxon read proportion allowed to retain a sample taxon
     marker : str 
          The sample marker name to search
     
     Returns
     -------
          list
               a list of pandas dataframe objects with marker taxon class read
               abundance proportions
     """

     print(f"{utils.get_localtime()} - - - - Computing {tax_class} taxonomy class read proportions for {marker} marker...\n")
     
     # load marker vsearch results and taxonomy lineage tables
     blast6out_names = ["query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", 
                         "evalue", "bits"]
     blast6out = pd.read_csv(align_resluts, names=blast6out_names, sep="\t")
     lineage_names = ["taxon_id", "order", "family", "genus", "species"]
     lineage = pd.read_table(marker_lineage, names=lineage_names, sep="\t") 
     merged_tables = pd.merge(blast6out, lineage.drop_duplicates(subset=['taxon_id'], keep='first'),
                         how="left", left_on="target", right_on="taxon_id")

     # select taxonomy class, quantify taxon read proportions,
     # and filter out taxonomic classes with low proportions
     merged_tables["reads"] = 1
     merged_tables = merged_tables[[tax_class, "reads"]]
     read_counts = merged_tables.groupby(tax_class).count().reset_index()
     read_counts[marker] = read_counts["reads"] / merged_tables["reads"].count()
     read_counts = read_counts[[tax_class, marker]]
     read_counts = read_counts[read_counts[marker] >= min_proportion]
     return read_counts

def get_rescaled_median_proportion(marker_dfs, tax_class, max_markers):
     """Computing and rescaling median proportions for the taxonomy class across all markers

     Parameters
     ----------
     marker_dfs : list
          a list of pandas dataframe objects with marker taxon class read
          abundance proportions
     tax_class : str
          Taxonomy class for quantify taxon level marker read abundance
     max_markers : int
          Maximum missing markers allowed to retain a sample taxon
     
     Returns
     -------
          pandas dataframe object
               a pandas dataframes object with marker taxon class read
               abundance proportions, median proportions, and rescaled
               median proportions
     """

     print(f"{utils.get_localtime()} - - - - Computing and rescaling median proportions for {tax_class} taxonomy class across all markers...\n")

     # compute rescaled median taxon reads proportions
     sample_df = reduce(lambda x, y: pd.merge(x, y, how="outer", on=[tax_class]), marker_dfs)
     max_markers += 1
     sample_df = sample_df[sample_df.isnull().sum(axis=1) < max_markers].fillna(0)
     sample_df["median"] = sample_df.iloc[:,1:].median(axis=1)
     sample_df["proportion"] = sample_df["median"] / sample_df["median"].sum()
     return sample_df

def main():
     t0 = time.time()
     args = utils.read_parameters("classify")

     print(f"{get_localtime()} - Starting read classification...\n")
     if args.min_proportion >= 0.01:
          raise Exception(f"The minimum taxon read proportion allowed to retain a sample taxon cannot exceed 0.1%!")
     search_markers(args.FASTA_DIR, args.DB_DIR, args.CONFIG_FILE, args.vsearch_aligner, args.tax_class, args.min_proportion, args.max_markers, args.threads)
 
     t1 = time.time()
     print(f"{get_localtime()} - Completed read classification...")
     print(f"Total elapsed time {int(t1 - t0)}\n")
     sys.exit(0)

if __name__ == '__main__':
    main()
