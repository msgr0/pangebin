import sys
import os
import pandas as pd
import numpy as np
import argparse as ap
import gfapy as gf


"""
Modifies output bins of Plasbin-Flow applied on a pangenome graph in order to account for contig fragmentation.
"""


def main(): 
    argparse = ap.ArgumentParser(description="Modify pangenome bins based on the input graph")
    argparse.add_argument("-p", "--pred", help="Input bins file", required=True)
    argparse.add_argument("-o", "--out", help="Output bins file", required=True)

    argparse.add_argument("--naive", help="Use naive super-sets", action="store_true")
    argparse.add_argument("--n", help="TYPE B EVAL", type=int)
    
    argparse.add_argument("--super", help="Use super-sets", action="store_true")
    argparse.add_argument("--graph", help="Input graph file")
    
    args = argparse.parse_args()




    bins = pd.read_csv(args.pred, sep="\t")
 
    bin_set = {}
    for i, row in bins.iterrows():
        if row["plasmid"] not in bin_set.keys():
            bin_set[row["plasmid"]] = []
            bin_set[row["plasmid"]].append(row["contig"])
        else:
            bin_set[row["plasmid"]].append(row["contig"])

    # NAIVE SUPER-SETS  EVAL(B.2)
    if (args.naive):
        new_bins = set()
        for bin in bin_set:
            candidates = []
            for bin_2 in bin_set:
                if bin_2 != bin:
                    if len(set(bin_set[bin]).intersection(set(bin_set[bin_2]))) > 0:
                        candidates.append(bin_2)
            new_bins.add(tuple(candidates))

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
            for contig in new_bin_set[bin]:
                bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(bin), "contig":  str(contig), "contig_len": bins.loc[bins["contig"] == contig, "contig_len"].values[0]}

        if (args.n == 1):
            bin_out_csv.to_csv(args.out, sep="\t", index=False)
            sys.exit(0)
        elif (args.n == 2):
            combined_bins = pd.concat([bins, bin_out_csv], ignore_index=True)
            combined_bins.to_csv(args.out, sep="\t", index=False)
            sys.exit(0)
    

    if (args.super):
        if (args.graph == None):
            print("Please provide a graph file")
            sys.exit(1)
        # graph in gfa format
        gfa = gf.Gfa.from_file(args.graph)

        ## SUPER-SETS EVAL(C)
        ## recover the set of fragments in each bin (bin_set)
        ## from the pangenome.gfa extract the set of fragment for each contig
        ## and store in a dictionary
        contig_set = {}
        for segment in gfa.segments:
            contig_set[segment.name] = []
            contigs = segment.cl.split(',')
            for contig in contigs:
                contig_set[segment.name].append(contig)

        # print(contig_set)

        fragment_set = {}
        for line in gfa.lines:
            if line.record_type == "P":
                fragment_set[line.name] = [str(x)[:-1] for x in line.segment_names]
        # print(fragment_set)

        ## 1. for each bin in bin_set, new_bins contains the set of contigs in contig_set[bin[fragment]]
        new_bins = set()
        for bin in bin_set:
            contigs = []
            for fragment in bin_set[bin]:
                contigs+=[x for x in contig_set[str(fragment)]]
            contigs = set(contigs)
            new_bins.add(tuple(set(contigs)))

        ordered_new_bins = sorted(new_bins, key=len)
        # for b in ordered_new_bins:
        #     print(len(b))
        print("ordered_new_bins", ordered_new_bins)
        print("")

        fragment_bins = []
        for b in ordered_new_bins:
            bin_contigs = set()
            for contig in b:
                for x in fragment_set[contig]:
                    bin_contigs.add(x)
            fragment_bins.append(bin_contigs)
        fragment_bins = sorted(fragment_bins, key=len) 
        print("fragment_bins", fragment_bins)
        print("")


        counter = 1
        fragment_bin_set = {}
        for bin in fragment_bins:
            fragment_bin_set[f"B{counter}"] = bin
            counter += 1
        print("fragment_bin_set", fragment_bin_set)
        print("")


        # build the bin_fragment_set, with each fragment as key and a bin set as value
        bin_fragment_set = {}
        for bin in fragment_bin_set:
            for fragment in fragment_bin_set[bin]:
                if fragment not in bin_fragment_set.keys():
                    bin_fragment_set[fragment] = set()
                    bin_fragment_set[fragment].add(bin)
                else:
                    bin_fragment_set[fragment].add(bin)
        print("bin_fragment_set", bin_fragment_set)
        print("")

        
        superbins_set = set()
        for x in bin_fragment_set:
            superbins_set.add(tuple(bin_fragment_set[x]))
        print("superbins_set", len(superbins_set),  superbins_set)
        print("")

        # remove bins that are subsets of other bins
        to_remove = set()
        for bin in superbins_set:
            for bin_2 in superbins_set:
                if bin != bin_2:
                    if set(bin).issubset(set(bin_2)):
                        to_remove.add(bin)
        superbins_set = superbins_set.difference(to_remove)
        print("superbins_set", len(superbins_set),  superbins_set)
        print("")

        # merge bins in superbins_set that have a non-empty intersection
        for bin in superbins_set.copy():
            to_remove = set()
            to_merge = set()
            for bin_2 in superbins_set.copy():
                if bin != bin_2:
                    if len(set(bin).intersection(set(bin_2))) > 0:
                        to_remove.add(bin)
                        to_remove.add(bin_2)
                        to_merge = to_merge.union(set(bin).union(set(bin_2)))
            superbins_set = superbins_set.difference(to_remove)
            if len(to_merge) > 0:
                superbins_set.add(tuple(to_merge))
        
        print("superbins_set", len(superbins_set),  superbins_set)
        print("")
        
        # csv output
        bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
        counter = 1
        for superbin in superbins_set:
            for bin in superbin:
                for contig in fragment_bin_set[bin]:
                    # print("contig", contig)
                    contig_len = len(gfa.segment(contig).sequence)
                    bin_out_csv.loc[len(bin_out_csv)] = {"plasmid": str(f"B{counter}"), "contig":  str(contig), "contig_len": contig_len }
            counter += 1
        pd.DataFrame(bin_out_csv).to_csv(args.out, sep="\t", index=False)



if __name__ == "__main__":
    main()
    
    
    

