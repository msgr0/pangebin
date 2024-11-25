from Bio import SeqIO as seq
from Bio import Seq
import sys
import argparse as ap
from itertools import combinations


"""
Remove FASTA entries based on a threshold sequence length.
"""




def main(args):
    fasta_input = args.input
    fasta_dict = dict()
    with open(fasta_input, "r") as spe:
        for record in seq.parse(spe, "fasta"):
            if len(record) >= args.threshold:
                fasta_dict[record.id] = record.seq
    with open(args.output, "w") as out:
        for x, rec in fasta_dict.items():
            out.write(f">{x}\n")
            out.write(f"{str(rec)}\n")


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="remove gfa segments below the threshold and connect neighbours"
    )
    parser.add_argument("-i", "--input", help="input graph to trim")
    parser.add_argument("-o", "--output", help="output graph")
    parser.add_argument(
        "-t", "--threshold", help="threshold to remove contigs", type=int
    )
    parser.add_argument("-c", "--contains", help="string to remove contigs")

    args = parser.parse_args()

    main(args)
