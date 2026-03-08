#!/usr/bin/env python3
"""
Filter for DMSP hmmsearch '--tblout' output.

Outputs a tab-delimited table with:
sample, gene_ID, dmsp_ID

Rules:
- For each gene across multiple hmmout files, keep only the best DMSP hit by:
    1) smaller e-value wins
    2) if e-value ties, larger score wins

Usage:
    python dmsp_filter.py -i <dir_with_hmmouts> -o <output.tsv> -s <sample_name>
"""

import argparse
import os
import glob
import math

def main():
    parser = argparse.ArgumentParser(description="Filter hmmsearch tblout results for DMSP")
    parser.add_argument("-i", "--input-dir", required=True, help="Directory containing .hmmout files")
    parser.add_argument("-o", "--output-file", required=True, help="Output filtered TSV")
    parser.add_argument("-s", "--sample", required=True, help="Sample name")
    
    args = parser.parse_args()
    
    best_hits = {} # gene_id -> (dmsp_id, evalue, score)
    
    hmmout_files = glob.glob(os.path.join(args.input_dir, "*.hmmout"))
    
    for hmmout_file in hmmout_files:
        if not os.path.exists(hmmout_file):
            continue
            
        with open(hmmout_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                    
                gene_id = parts[0]
                dmsp_id = parts[2]
                
                try:
                    evalue = float(parts[4])
                    score = float(parts[5])
                except ValueError:
                    continue
                
                # Update best hit
                if gene_id not in best_hits:
                    best_hits[gene_id] = (dmsp_id, evalue, score)
                else:
                    curr_evalue = best_hits[gene_id][1]
                    curr_score = best_hits[gene_id][2]
                    
                    if evalue < curr_evalue:
                        best_hits[gene_id] = (dmsp_id, evalue, score)
                    elif evalue == curr_evalue and score > curr_score:
                        best_hits[gene_id] = (dmsp_id, evalue, score)
                        
    # Write output
    with open(args.output_file, "w", encoding="utf-8") as out:
        for gene_id in sorted(best_hits.keys()):
            dmsp_id = best_hits[gene_id][0]
            out.write(f"{args.sample}\t{gene_id}\t{dmsp_id}\n")

if __name__ == "__main__":
    main()
