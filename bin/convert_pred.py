import os
import numpy as np
import pandas as pd
import argparse
import re


def main():
    parser = argparse.ArgumentParser(description='Convert prediction to plaseval')
    parser.add_argument('--input', help='path to prediction mixed file')
    parser.add_argument('--output-gt', help='path to submission file')
    parser.add_argument('--output-bin', help='path to submission file')
    args = parser.parse_args()

    out_bin_path = args.output_bin
    out_gt_path = args.output_gt
    
    df = pd.read_csv(args.input)
    header = ['plasmid', 'contig','contig_len']
    out_gt = pd.DataFrame(columns=header)
    for index, row in df.iterrows():
        # Access the values of each column in the row
        plasmids = []
        if row['plasmids'] is not np.nan and row['chromosomes'] is np.nan:
            plasmids = row['plasmids'].split(',')
        contig = row['contig']
        contig_len = row['length']
        
        for plasmid in plasmids:
            out_gt.loc[len(out_gt)] = [plasmid, contig, contig_len]
            
    out_gt.to_csv(out_gt_path, index=False, sep='\t')
    
    out_bin = pd.DataFrame(columns=header)
    for index, row in df.iterrows():
        bins = []
        if row['bins'] is not np.nan:
            print('row_bins:',row['bins'])
            bins = ["P"+x for x in re.split('P', str(row['bins']))[1:]]
            print("bins:", bins)
        contig = row['contig']
        contig_len = row['length']
        for bin in bins:
            out_bin.loc[len(out_bin)] = [bin, contig, contig_len]
    
    out_bin.to_csv(out_bin_path, index=False, sep='\t')
    print('okj')
    
if __name__ == "__main__":
    main()