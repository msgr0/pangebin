import pandas as pd
import numpy as np
import os
import sys
import gfapy as gfa
import argparse
# import fastaparser as fp


def main():
    parser = argparse.ArgumentParser(description='Transform Gene Density file into GPLAS prediction')
    parser.add_argument('--input', help='Input file')
    parser.add_argument('--gfa', help='GFA file')
    
    parser.add_argument('--output', help='Output file')
    
    
    
    args = parser.parse_args()

    gene_density = pd.read_csv(args.input, sep='\t', header=None, skip_blank_lines=True)
    
    prediction_header = ['Prob_Chromosome', 'Prob_Plasmid', 'Prediction', 'Contig_name', 'Contig_length']
    prediction = pd.DataFrame(columns=prediction_header)
    gfa_file = gfa.Gfa.from_file(args.gfa)
    
    # seg_ids = []
    # with open(args.fasta) as fasta_file:
    #     parser = fp.Reader(fasta_file)
    #     for seq in parser:
    #         seg_ids.append(seq.id)
    
    for _, row in gene_density.iterrows():
        contig_name = row[0]
        plasmid_score = row[1]
        contig_length = gfa_file.segment(str(contig_name)).LN if gfa_file.segment(str(contig_name)).LN != None else len(gfa_file.segment(str(contig_name)).sequence)
        prob_plasmid = float(plasmid_score)
        prob_chromosome = 1.0 - prob_plasmid
        label = 'Plasmid' if prob_plasmid > prob_chromosome else 'Chromosome'
        # contig_name = seg_ids[index]
        prediction.loc[len(prediction)] = {'Prob_Chromosome': prob_chromosome, 'Prob_Plasmid': prob_plasmid, 'Prediction': label, 'Contig_name': contig_name, 'Contig_length': contig_length}
        
    prediction.to_csv(args.output, sep='\t', index=False)
    
if __name__ == "__main__":
    main()