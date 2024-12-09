import os
import sys
import numpy as np
import pandas as pd
import argparse as ap
import gfapy as gfa


def main(args):
    df_header = ["number", "Bin"]
    df = pd.read_csv(args.input, sep=" ")
    if args.prefix != None:
        assert args.prefix in ["uni", "ske", ""]
    else:
        args.prefix = ""

    header = ["plasmid", "contig", "contig_len"]
    out_prediction = pd.DataFrame(columns=header)

    graph = gfa.Gfa.from_file(args.gfa)

    for index, row in df.iterrows():
        bin_name = str(row["Bin"])
        if bin_name.startswith("U"):  # --> unbinned contig
            plasmid = "U"
        elif bin_name.startswith("I"):  # ---> isolated
            plasmid = "IS" + str(row["Bin"]).split("_")[1]
        else:
            try:
                plasmid = "P" + bin_name
            except:
                raise Exception("Could not retrive the bin name-type")
                sys.exit(1)
        contig = str(row["number"])
        out_contig = str(args.prefix) + contig
        contig_len = -1
        try:
            contig_len = graph.line(contig).LN
            if contig_len <= 0:
                raise Exception("contig.LN == 0")
        except:  ## contig len == 0 or other error catched in e1
            try:
                contig_len = len(graph.line(contig).sequence)
                if contig_len <= 0:
                    raise Exception("seq(contig).length == 0")
                elif contig_len > 0:
                    # print("warning: ")
                    pass
            except:
                print("Halted cause of length error: check your gfa file.")
                print("contig", contig)
                sys.exit(1)

        out_prediction.loc[len(out_prediction)] = [plasmid, out_contig, contig_len]

    out_prediction.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="Convert PLASBIn-FLOW prediction to PlasEval."
    )

    parser.add_argument("--input", "-i", help="Input prediction file of PBF")
    parser.add_argument("--gfa", "-g", help="Input graph that led to the prediction")
    parser.add_argument("--output", "-o", help="Output file, plaseval prediction")
    parser.add_argument(
        "--prefix",
        "-p",
        help="prefix to renaming the contigs, use with 'ske' or 'uni' for convert assebmly prediction",
    )

    main(parser.parse_args())
