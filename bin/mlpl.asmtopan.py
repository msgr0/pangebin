import argparse
import gfapy as gp
import numpy
import pandas as pd
import sys


"""
Converts prediction of mlplasmid for a pangenome graph
"""



def main(args):
    pred = args.pred
    graph = args.graph
    output = args.output
    pbf = args.pbf
    thr = args.thr if args.thr != None else 0
    print("starting converting MLPLSAMIDS prediction ...")

    mlpred_head = [
        "Prob_Chromosome",
        "Prob_Plasmid",
        "Prediction",
        "Contig_name",
        "Contig_length",
    ]

    mlpred = pd.read_csv(pred, sep="\t", header=0)
    mlpred["Contig_name"] = mlpred["Contig_name"].astype(str)

    # return
    df = pd.DataFrame(columns=mlpred_head)

    df_pbf = pd.DataFrame(columns=["contig", "pscore"])

    gfa = gp.Gfa.from_file(graph)

    for frag in gfa.segments:
        contig_list = frag.cl.split(",")  ## ogni frammento ha per forza una contig list
        ncontigs = len(contig_list)
        prob_chromo = 0
        prob_plas = 0
        length = frag.LN if frag.LN != None else len(frag.sequence)
        cumulated_l = 0
        for contig in contig_list:
            # print(contig)
            try:

                mlpred.loc[mlpred["Contig_name"] == contig, "Prob_Chromosome"].values[0]

            except:
                ncontigs -= 1
                continue

            prob_chromo += (
                mlpred.loc[mlpred["Contig_name"] == contig, "Prob_Chromosome"].values[0]
                * mlpred.loc[mlpred["Contig_name"] == contig, "Contig_length"].values[0]
            )
            prob_plas += (
                mlpred.loc[mlpred["Contig_name"] == contig, "Prob_Plasmid"].values[0]
                * mlpred.loc[mlpred["Contig_name"] == contig, "Contig_length"].values[0]
            )
            # labeling = mlpred.loc[mlpred[3] == contig, 2].values[0]
            cumulated_l += mlpred.loc[
                mlpred["Contig_name"] == str(contig), "Contig_length"
            ].values[0]
        if ncontigs == 0:
            continue

        prob_chromo /= cumulated_l
        prob_plas /= cumulated_l

        if prob_plas >= prob_chromo:
            label = "Plasmid"
        else:
            label = "Chromosome"

        name = f"S{frag.name}_LN:i:{frag.LN}_dp:f:{frag.cv}"

        if length < thr:
            continue

        if output != None:
            df.loc[len(df)] = {
                "Prob_Chromosome": prob_chromo,
                "Prob_Plasmid": prob_plas,
                "Prediction": label,
                "Contig_name": name,
                "Contig_length": length,
            }

        if pbf != None:
            df_pbf.loc[len(df_pbf)] = {"contig": frag.name, "pscore": prob_plas}

    if pbf != None:
        df_pbf.to_csv(pbf, sep="\t", index=False, header=False)
    if output != None:
        df.to_csv(output, sep="\t", index=False)
    print("done transforming ML PLASMID PRED")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert mlplasmid pred")

    parser.add_argument("--pred", "-p", help="mlplasmid asm pred")
    parser.add_argument("--graph", "-g", help="pangenome_graph convert mlplas into")
    parser.add_argument("--output", "-o", help="output")
    parser.add_argument("--pbf", help="pbf output mode")
    parser.add_argument("--thr", help="threshold", type=int)
    args = parser.parse_args()
    main(args)
