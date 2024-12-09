import argparse as ap
import re
import sys
import os
import fileinput

"""
Remove Skesa Assembly Graph branching endings based
on the skesa-provided naming conventions
"""


def main(input, output):
    outfile = open(file=output, mode="w")
    with open(file=input, mode="r") as infile:
        for line in infile:
            if line.startswith("S"):
                if re.search("^S.*:\d.:.*\t.*", line) == None:
                    outfile.write(
                        line
                    )  # this is a skesa branching line we want to remove from our file
            if line.startswith("L"):
                outfile.write(re.sub(":\d.:\d*:\d*\t", "\t", line))
    outfile.close()

    lines = open(output, "r").readlines()
    lines_set = set(lines)
    out = open(output, "w")
    for line in lines_set:
        out.write(line)
    out.close()


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input", help="path to the graph to trim")
    parser.add_argument("-o", "--output", help="path where to store the outfile")
    args = parser.parse_args()
    main(args.input, args.output)
