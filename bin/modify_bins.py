import sys
import os
import pandas as pd
import altair as alt
import numpy as np
import argparse as ap
import gfapy as gfa

"""
Modifies output bins of Plasbin-Flow applied on a pangenome graph in order to account for contig fragmentation.
"""

def main(): 
    argparse = ap.ArgumentParser(description="Modify pangenome bins based on the input graph")
    argparse.add_argument("-g", "--graph", help="Input graph file")
    argparse.add_argument("-p", "--gtruth", help="Input gtruth file", required=True)
    argparse.add_argument("-b", "--bins", help="Input bins file", required=True)
    argparse.add_argument("-o", "--output", help="Output bins file", required=True)
    argparse.add_argument("-m", "--mod", help="Mod value for bin modification", required=True, type=int)
    args = argparse.parse_args()
    
    # Load graph
    # graph = gfa.Gfa.from_file(args.graph)
    # Load bins
    bins = pd.read_csv(args.bins, sep="\t")
    # Modify bins according to mod 1
    truth = pd.read_csv(args.gtruth, sep="\t")
    
    ### get bins ###
    bin_set = {} # {B1: [contig1, contig2, ...], B2: [contig3, contig4, ...]}
    for i, row in bins.iterrows():
        if row["plasmid"] not in bin_set.keys():
            bin_set[row["plasmid"]] = []
            bin_set[row["plasmid"]].append(row["contig"])
        else:
            bin_set[row["plasmid"]].append(row["contig"])
    

    
    if args.mod == 1:
        new_bins = set()
        for bin in bin_set:
            candidates = []
            for bin_2 in bin_set:
                if len(set(bin_set[bin]).intersection(set(bin_set[bin_2]))) > 0:
                    candidates.append(bin_2)
            new_bins.add(tuple(candidates))
            
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
            
        bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
        for bin in new_bin_set:
            bin_len = 0
            for contig in new_bin_set[bin]:
                contig_len = truth.loc[truth["contig"] == contig, "contig_len"]
                if len(contig_len) > 0:
                    contig_len = contig_len.values[0]
                else:
                    contig_len = -1                    
                bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": contig_len}
        bin_out_csv.to_csv(args.output, sep="\t", index=False)
    
    
    elif args.mod == 2:
        graph = gfa.Gfa.from_file(args.graph)
        
        # augment each bin with the fragments of the pangenome that belongs to the same contig
        new_bins = set()
        for bin in bin_set:

            contigs = set() ## the fragments of contigs to add to the bin
            fragments = set(bin_set[bin])
            for f in fragments:
                seg = graph.segment(f)
                if (seg != None):
                    contig_list = seg.cl.split(",")
                    for c in contig_list:
                        contigs.add(c)
                        
            candidates = set()
            for c in contigs:
                path = graph.line(str(c))
                if path != None:
                    segs = path.segment_names
                    for seg in segs:
                        candidates.add(str(seg.name))
            new_bins.add(tuple(candidates))
                        
        new_bin_set = {}
        counter = 1
        for bin in new_bins:
            contigs = list(set(bin))
            new_bin_set[f"B{counter}"] = contigs
            counter += 1
            
            
        bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
        for bin in new_bin_set:
            bin_len = 0
            for contig in new_bin_set[bin]:
                contig_len = truth.loc[truth["contig"] == contig, "contig_len"]
                if len(contig_len) > 0:
                    contig_len = contig_len.values[0]
                else:
                    contig_len = -1                    
                bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": contig_len}
        bin_out_csv.to_csv(args.output, sep="\t", index=False)
    
    else:
        pass
        # mod_1()
        # mod_2()
    

if __name__ == "__main__":
    main()
    
    
    

