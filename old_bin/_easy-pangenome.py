import argparse
import gfapy as gp
import re
import os
import sys


"""
Pangenome to GPLAS pangenome. Remove Paths and convert the pangenome to 'unicycler' format.
"""

def main(input, output):
    print("Input: ", input)
    print("Output: ", output)

    gfa = gp.Gfa.from_file(input)
    gfa.validate()
    for line in gfa.lines:
        if str(line)[0] == "P":
            line.disconnect()
        if str(line)[0] == "L":
            line.aa = None
            line.lt = None
        if str(line)[0] == "S":
            # new_name = "S" + line.name + "__"
            # line.name = new_name
            line.OC = None
            line.cl = None
            line.aa = None
            line.ap = None
            line.dp = line.cv
            line.cv = None
    gfa.to_file(f"{str(output)}")

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pangenome to GPLAS pangenome")
    parser.add_argument("--input", help="Input gfa pangenome file")
    parser.add_argument("--output", help="Output file")

    args = parser.parse_args()
    if args.output == None:
        main(args.input, f"{str(args.input)[:-4]}.slim.gfa")
    else:
        main(args.input, args.output)
