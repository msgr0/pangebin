import gfapy as gp
import argparse as ap

"""
Rename GFA contigs with cardinal numbers.
Convert also KC of skesa assembly graphs to dp of uncycler graphs
"""

def main(args):
    graph = gp.Gfa.from_file(args.input, vlevel=0)
    graph.validate()

    if args.prefix != None:
        rename_gfa(graph, args.prefix)
        graph.validate()

    if args.convert:
        convert(graph)
        graph.validate()

    graph.to_file(args.output)


def convert(gfa):
    total_coverage = 0
    total_length = 0
    for seg in gfa.segments:
        total_coverage += seg.KC
        if seg.LN == None:
            seg.set_datatype("LN", "i")
            seg.LN = len(seg.sequence)
        total_length += seg.LN

    for seg in gfa.segments:
        seg.set_datatype("dp", "f")
        seg.dp = float((seg.KC * total_length) / (seg.LN * total_coverage))
        seg.KC = None


def rename_gfa(gfa, prefix):
    gfa.vlevel = 3
    counter = 1
    for seg in gfa.segments:
        if seg.LN is not None:
            if seg.LN == 0:
                seg.sequence = "*"
                seg.LN = 0
                gfa.validate()
        else:
            if len(seg.sequence) == 0:
                seg.sequence = "*"
                seg.LN = 0
                gfa.validate()
        seg.name = f"{prefix}{counter}"
        counter += 1


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Rename a GFA file with 'prefix'")
    parser.add_argument("-i", "--input", help="Input GFA file")
    parser.add_argument("-o", "--output", help="Output GFA file")
    parser.add_argument("-p", "--prefix", help="Type of the GFA file")
    parser.add_argument(
        "-c",
        "--convert",
        help="Convert coverage if skesa or pangenome into unicycler",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)
