## build a CSV ground truth using blast of a PAN_ASSEMBLY graph

## FNA of the genome we want to evaluate
## INPUT of the pangenome we want to label
## INPUT of the REFERENCE we want to build the GT on (unicycler or skesa graph)

## OUTPUT reference labeling, to evaluate directly the assembly
## OUTPUT input labeling, to evaluate the pangenome graph in respect of the assembly

## IF MIXED, PROVIDE ALL THE INPUT ASSEMBLY AND RUN IN THE ORIGINAL MODE.


# steps:
# 1.
import os
import sys
import numpy as np
import pandas as pd
import argparse as ap
import gfapy as gpy

from mappings_utils import run_blast6, DEFAULT_PID_THRESHOLD, DEFAULT_COV_THRESHOLD

from log_errors_utils import (
    check_file,
    log_file,
    CustomException,
    clean_files,
    create_directory,
    run_cmd,
    process_exception,
)
from ground_truth import compute_ground_truth_file
from gfa_fasta_utils import read_FASTA_len, read_FASTA_id

pid_threshold = DEFAULT_PID_THRESHOLD
cov_threshold = DEFAULT_COV_THRESHOLD


def convert_gt(pangenome_graph, gt_csv_file, pangenome_gt_file, mix):

    header_def = ["pls_id", "ctg_id", "coverage_str", "pls_len", "ctg_len"]
    gt_csv = pd.read_csv(gt_csv_file, sep="\t", header=None)
    out_csv = pd.DataFrame(columns=header_def)

    pangenome = gpy.Gfa.from_file(pangenome_graph)
    pangenome.validate()

    prefix = (gt_csv[1].values[0])[0:2]

    for seg in pangenome.segments:
        contig_list = seg.cl.split(",")
        plasmid_label = None
        label = None
        coverage = 0
        length = 0

        for contig in contig_list:
            if not contig.startswith(prefix) and not mix:
                continue

            if (
                contig in gt_csv[1].values
            ):  # if the contig is in the ground truth file, OVVERO if the contig has been matched against a PLASMID
                new_label = gt_csv.loc[gt_csv[1] == contig, 0].values[0]
                if new_label.startswith("pl"):
                    plasmid_label = (
                        "mixed"
                        if plasmid_label != None and plasmid_label != new_label
                        else new_label
                    )
                    new_label = "plasmid"
                elif new_label.startswith("ch"):
                    new_label = "chromosome"
                else:
                    assert False, f"Should be either plasmid or chromosome. Check FNA"
                if length == 0:
                    length = gt_csv.loc[gt_csv[1] == contig, 3].values[0]
                coverage = max(coverage, gt_csv.loc[gt_csv[1] == contig, 2].values[0])
                if label == "chromosome" or label == None:
                    label = new_label
                elif label == "plasmid":
                    pass

                # elif label == "chromosome":
                #     pass
                # elif label == "plasmid":
                #     pass

        """Add the new segment(fragment) information to the transformed ground truth"""
        # header = ["pls_id", "ctg_id", "coverage_str", "pls_len", "ctg_len"]
        if label == "plasmid":
            out_csv.loc[len(out_csv)] = {
                "pls_id": plasmid_label,
                "ctg_id": seg.name,
                "coverage_str": coverage,
                "pls_len": length,
                "ctg_len": seg.LN,
            }
        elif label == "chromosome":
            out_csv.loc[len(out_csv)] = {
                "pls_id": "chromosome",
                "ctg_id": seg.name,
                "coverage_str": coverage,
                "pls_len": length,
                "ctg_len": seg.LN,
            }
        elif label == None:
            print("skipping fragment \t", seg.name, "\tits not related to asmb")
        out_csv.to_csv(pangenome_gt_file, sep="\t", header=None, index=None)

    # for index, row in gt_csv.iterrows():
    #     # row[1] va cercato in gfa.Paths; per ogni fragment dentro il path, appendo
    #     for frag in frag_list:
    #         out_csv.loc[len(out_csv)] = {row[0], frag.name, row[2], row[3], row[4]}

    #             if contig in ground_truth['contig'].values: # if the contig is in the ground truth file
    #                 new_label = ground_truth.loc[ground_truth['contig'] == contig, 'label'].values[0] # retrieve the label
    #                 if new_label in ["plasmid", "chromosome", "ambiguous"]:
    #                     if label is None or label == "undef": # if this is the first label, assign it
    #                         label = new_label
    #                     elif label != new_label: # if it is not the first label, is ambiguous if it differs from the previous one
    #                         label = "ambiguous"
    #                     elif label == new_label: # (this assignment is not necessary, but it is more explicit)
    #                         label = new_label
    #                 else:
    #                     continue # do not annotate [short] and [undef] label
    #             else:
    #                 assert(False), f"Contig {contig} not found in the ground truth file"

    #             if new_label is None:
    #                 label = "undef"

    #         """Add the new segment(fragment) information to the transformed ground truth"""
    #         gt_transformed.loc[len(gt_transformed)] = {'contig': seg.name, 'label': label, 'length': seg.LN, 'type': seg.aa}


def main(args):
    # args.pangenome
    # args.assembly
    # args.reference
    mapping_file = args.output + ".mapping.tsv"
    assembly_gt = args.output + ".gt.tsv"
    pangenome_gt = args.output + ".pan.gt.tsv"
    

    # split fasta based on provenience, plasmid-chromosome?? or just strip the chromosome

    # The following computes the ASSEMBLY GT
    run_blast6(args.assembly, args.reference, mapping_file)
    compute_ground_truth_file(
        args.output,
        mapping_file,
        read_FASTA_len(args.assembly, gzipped=False),
        read_FASTA_len(args.reference, gzipped=False),
        0.95,
        0.5,
        assembly_gt,
    )

    mix = args.mix if args.mix else False
    convert_gt(args.pangenome, assembly_gt, pangenome_gt, mix)


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="Convert PLASBIn-FLOW prediction to PlasEval."
    )

    parser.add_argument(
        "--pangenome", "-p", help="Input pangenome graph built from Skesa and Unicycler"
    )
    parser.add_argument(
        "--mix",
        "-m",
        help="if mixed mode, provide the plasmid-mixed pangenome graph (MIXED)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--assembly",
        "-a",
        help="Input assembly to base the gt into (skesa, unicycler or MIXED)",
    )
    parser.add_argument(
        "--reference", "-b", help="Input genomic.fna retrieved from NCBI"
    )
    parser.add_argument(
        "--output",
        "-o",
        help="SAMPLE_ID: will be used as the output PREFIX, for 1.mapping plasmids file, 2.pangenome GT and 3.ASSEMBLY GT.",
    )
    args = parser.parse_args()
    if args != None:
        main(args)
