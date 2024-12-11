from dataclasses import dataclass
from typing import Annotated
from pathlib import Path

import gfapy as gp # type: ignore[import]
import typer

APP = typer.Typer(rich_markup_mode="rich")

# XXX: Add Gplas Utils to extend pangebin with Gplas solver
# @dataclass
# class GplasUtilsArgs:

# TODO add pbf utils to convert bins to single column prediction
# TODO add to pbf utils to convert bins to ovl bins and/or naive bins
# @dataclass
# class PbfUtilsArgs:


@dataclass
class UtilsArgs:
    """GFAUtils arguments"""

    ARG_INPUT_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_INPIUT_FASTA = typer.Argument(
         help="Input FASTA file",
    )

    ARG_OUTPUT = typer.Argument(
        help="Output file name",
    )

    OPT_PREFIX = typer.Option(
        "--prefix",
        help="Type of the GFA file",
    )

@APP.command()
def rename_contigs(
    gfa: Annotated[Path, UtilsArgs.ARG_INPUT_GFA],
    prefix: Annotated[str, UtilsArgs.OPT_PREFIX],
)-> None:

    graph = gp.Gfa.from_file(gfa, vlevel=0)
    graph.validate()
    graph.vlevel = 3
    counter = 1
    for seg in graph.segments:
        if seg.LN is not None:
            if seg.LN == 0:
                seg.sequence = "*"
                seg.LN = 0
                graph.validate()
        else:
            if len(seg.sequence) == 0:
                seg.sequence = "*"
                seg.LN = 0
                graph.validate()
        seg.name = f"{prefix}{counter}"
        counter += 1
    graph.validate()
    graph.to_file(args.output)

@APP.command()
def convert_KC_to_dp(
    gfa: Annotated[Path, UtilsArgs.ARG_INPUT_GFA],
    output: Annotated[str, UtilsArgs.ARG_OUTPUT],
)-> None:
    graph = gp.Gfa.from_file(gfa, vlevel=0)
    graph.validate()
    total_coverage = 0
    total_length = 0
    for seg in graph.segments:
        total_coverage += seg.KC
        if seg.LN == None:
            seg.set_datatype("LN", "i")
            seg.LN = len(seg.sequence)
        total_length += seg.LN

    for seg in graph.segments:
        seg.set_datatype("dp", "f")
        seg.dp = float((seg.KC * total_length) / (seg.LN * total_coverage))
        seg.KC = None
    graph.validate()
    graph.to_file(output)




@APP.command()
def remove_nodes():
    pass


"""
REMOVE nodes remove_nodes.py
"""


@APP.command()
def extract_fasta():
    pass


"""
extract FASTA
"""


@APP.command()
def gfa_cleaner():
    pass


"""
(PANGENOME) GFA CLEANER
"""



@APP.command()
"""
remove contig if contains name (string_fasta_plasmid.py) """
def remove_contig():
    pass

@APP.command()
""" extract a gfa into fasta file """
def extract_gfa():
    pass



if __name__ == "__main__":
    APP()
