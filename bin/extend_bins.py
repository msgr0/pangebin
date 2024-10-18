import sys
import os
import pandas as pd
import numpy as np
import argparse as ap


def main(): 
    argparse = ap.ArgumentParser(description="Modify pangenome bins based on the input graph")
    argparse.add_argument("-p", "--pred", help="Input bins file", required=True)
    argparse.add_argument("-o", "--out", help="Output bins file", required=True)
    args = argparse.parse_args()

    bins = pd.read_csv(args.pred, sep="\t")
 
    bin_set = {}
    for i, row in bins.iterrows():
        if row["plasmid"] not in bin_set.keys():
            bin_set[row["plasmid"]] = []
            bin_set[row["plasmid"]].append(row["contig"])
        else:
            bin_set[row["plasmid"]].append(row["contig"])
    # print(bins)
    # print(bin_set)

    new_bins = set()
    for bin in bin_set:
        candidates = []
        for bin_2 in bin_set:
            if bin_2 != bin:
                if len(set(bin_set[bin]).intersection(set(bin_set[bin_2]))) > 0:
                    candidates.append(bin_2)
        new_bins.add(tuple(candidates))
    # print(new_bins)
    #     new_bins = set()
    #     for bin in bin_set:
    #         candidates = []
    #         for bin_2 in bin_set:
    #             if len(set(bin_set[bin]).intersection(set(bin_set[bin_2]))) > 0:
    #                 candidates.append(bin_2)
    #         new_bins.add(tuple(candidates))
            
    # mod_1()
    new_bin_set = {}
    counter = 1
    for bin in new_bins:
        contigs = []
        for i in bin:
            contigs+=[x for x in bin_set[i]]
        contigs = list(set(contigs))
        new_bin_set[f"B{counter}"] = contigs
        counter += 1

    # print(new_bin_set)
    bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
    for bin in new_bin_set:
        for contig in new_bin_set[bin]:
            bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": bins.loc[bins["contig"] == contig, "contig_len"].values[0]}

    # print(bin_out_csv)

    # concat bin_out_csv with bins
    combined_bins = pd.concat([bins, bin_out_csv], ignore_index=True)

    combined_bins.to_csv(args.out, sep="\t", index=False)


    # bin_out_csv.to_csv(args.out, sep="\t", index=False)
    #     bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
    #     for bin in new_bin_set:
    #         bin_len = 0
    #         for contig in new_bin_set[bin]:
    #             contig_len = truth.loc[truth["contig"] == contig, "contig_len"]
    #             if len(contig_len) > 0:
    #                 contig_len = contig_len.values[0]
    #             else:
    #                 contig_len = -1                    
    #             bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": contig_len}
    #     bin_out_csv.to_csv(args.output, sep="\t", index=False)
    
            
    #     bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
    #     for bin in new_bin_set:
    #         bin_len = 0
    #         for contig in new_bin_set[bin]:
    #             contig_len = truth.loc[truth["contig"] == contig, "contig_len"]
    #             if len(contig_len) > 0:
    #                 contig_len = contig_len.values[0]
    #             else:
    #                 contig_len = -1                    
    #             bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": contig_len}
    #     bin_out_csv.to_csv(args.output, sep="\t", index=False)
    
    # else:
    #     pass


if __name__ == "__main__":
    main()
    
    
    

