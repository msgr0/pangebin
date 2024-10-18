import os
import sys
import numpy
import pandas as pd
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(args):
    fasta_p = args.input
    fasta_o = args.output

    #  take every header, if contains no "plasmid" remove and print headers for debug
    fasta_open = open(fasta_p, "r")
    id_fun = lambda x: x
    record_fun = lambda x: x
    ctgs_dict = {}
    ctgs_seq = {}

    for record in SeqIO.to_dict(SeqIO.parse(fasta_p, "fasta")).items():
        id = record[0]
        descp = record[1].description[len(record[0]) + 1 :]
        seq = record[1].seq
        ctgs_dict.update({id: descp})
        ctgs_seq.update({id: seq})

    ctgs_chr = ctgs_seq.copy()
    for key in ctgs_dict:
        desc = str(ctgs_dict[key])
        if "plasmid" in desc:
            print("adding...", key, desc)
            ctgs_chr.pop(key)
            continue
        else:
            print("removing... ", key, desc)
            ctgs_seq.pop(key)
    with open(fasta_o, "w") as f:
        counter = 0
        for key in ctgs_chr:
            counter += 1
            f.write(">chromosome" + str(counter) + "\n")
            f.write(str(ctgs_chr[key]) + "\n")
        counter = 0
        for key in ctgs_seq:
            counter += 1
            # out = ">plasmid" + str(counter) + "\n"
            f.write(">plasmid" + str(counter) + "\n")
            f.write(str(ctgs_seq[key]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="strip fasta file to plasmids")
    parser.add_argument("--input", "-i", help="refseq genomic fna from ncbi")
    parser.add_argument("--output", "-o", help="output")

    args = parser.parse_args()
    main(args)
