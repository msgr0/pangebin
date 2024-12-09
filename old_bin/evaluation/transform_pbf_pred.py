import os
import sys
import numpy as np
import pandas as pd
import argparse as ap
import gfapy as gfa


def main(args):

    df_header = ["#Pls_ID", "Flow", "GC_bin", "Contigs"]
    df = pd.read_csv(args.input, sep="\t")

    header = ["plasmid", "contig", "contig_len"]
    out_prediction = pd.DataFrame(columns=header)

    graph = gfa.Gfa.from_file(args.gfa)

    for index, row in df.iterrows():
        plasmid = str(row["#Pls_ID"])
        contigs = []
        if (str(row["GC_bin"]) == "0"):
            continue
        for contig in str(row["Contigs"]).split(","):
            contigs.append(contig.split(":")[0])
        for contig in contigs:
            contig_len = -1
            try:
                contig_len = graph.line(contig).LN
                if contig_len <= 0:
                    raise Exception("contig.LN == 0")
            except Exception as e1:  ## contig len == 0 or other error catched in e1
                try:
                    contig_len = len(graph.line(contig).sequence)
                    if contig_len <= 0:
                        raise Exception("seq(contig).length == 0")
                    elif contig_len > 0:
                        print("warning: ", e1)
                except e2:
                    print(
                        "Halted cause of length error:", e1, e2, "check your gfa file."
                    )
                    sys.exit(1)
            out_prediction.loc[len(out_prediction)] = [plasmid, contig, contig_len]
    out_prediction.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="Convert PLASBIn-FLOW prediction to PlasEval."
    )

    parser.add_argument("--input", "-i", help="Input prediction file of PBF")
    parser.add_argument("--gfa", "-g", help="Input graph that led to the prediction")
    parser.add_argument("--output", "-o", help="Output file, plaseval prediction")

    main(parser.parse_args())
