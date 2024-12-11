from dataclasses import dataclass
from typing import Annotated
from pathlib import Path

import gfapy as gp # type: ignore[import]
import typer

APP = typer.Typer(rich_markup_mode="rich")

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


o
@APP.command()
def rename_contigs(
    gfa: Annotated[Path, UtilsArgs.ARG_INPUT_GFA],
    prefix: Annotated[str, UtilsArgs.OPT_PREFIX],
):

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


    # if args.convert:
    #     convert(graph)
    #     graph.validate()

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
