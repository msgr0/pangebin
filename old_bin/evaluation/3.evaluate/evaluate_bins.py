"""Methods:
- precise description of the protocole defining the ground truth, with the thresholds used
- for each plasmid binning tool used, same than above

Data/results, for each sample
- assembly graph:
 contribution (in actual cumulated length) of fragment of the pangenome graph that are (1) common to both assemblers, (2) specific to unicycler, (3) specific to Skesa; this can be shown as a single bar in a bar graph
- ground truth:
 - same but for all plasmids considered together
 - total length of the true plasmids
- binning results:
 - precision, recall, F1
 - for precision and recall, as above,  of fragment of the pangenome graph that are (1) common to both assemblers, (2) specific to unicycler, (3) specific to Skesa

cedric
"""

import argparse as ap

import pandas as pd

import altair as alt

# import gfapy as gp
import os
import numpy as np

# import . as utl


def safe_div(x, y):
    """Safe division to handle division by zero"""
    if y == 0:
        return 0
    else:
        return x / y


def safe_div(x, y):
    """Safe division to handle division by zero"""
    if y == 0:
        return 0
    else:
        return x / y


def parser():

    parser = ap.ArgumentParser()
    parser.add_argument(
        "--bin", help="Input bin file, tsv output of plasbin flow", required=True
    )
    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument(
        "--csv", help="Input csv file, ground truth file", required=True
    )
    # parser.add_argument("--gfa", help="Input gfa file containing the pangenome")
    parser.add_argument("--output", help="Output folder")
    parser.add_argument("--plot", help="Plot the results", action="store_true")

    return parser.parse_args()


def plot_(data: pd.DataFrame, _output: str, _title: str):
    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x="sample",
            y="length",
            color="type",
            order=alt.Order("type", sort="ascending"),
            tooltip="length",
        )
        .properties(title=_title)
    )
    chart.save(_output)


def main():
    """main function:
    - read the ground truth file
    - read the bin file, output of binning_tool, standardized [plasmid, contig, contig_len]
    - read the pangenome gfa file
    - calculate the statistics
    - print the results
    - plot the results
    """

    args = parser()

    sample = args.sample
    out = f"Sample: {sample}\n"
    """GROUND TRUTH FILE [header]
    plasmid_label, contig, coverage, plasmid_length, contig_length"""

    header_gt = [
        "plasmid_label",
        "contig",
        "coverage",
        "plasmid_length",
        "contig_length",
    ]
    ground_truth = pd.read_csv(args.csv, sep="\t", header=None)
    # gt_plasmid_length = pick one length, per type of plasmid
    # map {plasmid_name: length}
    print("ground truth actual plasmid length")

    # Read in the bin file
    header_bin = ["label", "contig", "contig_len"]
    bin_file = pd.read_csv(args.bin, sep="\t")

    # predicted_plasmid_length = sum(all the length in file gt, divided per plasmid_label)

    # bin_file.sort_values(by="contig", inplace=True)
    # bin_file = bin_file.groupby("contig").agg("sum").reset_index()
    # bin_file.drop_duplicates(subset="contig", keep="first", inplace=True)

    # pred_plasmids_len = bin_file.groupby("").sum()["length"]["plasmid"]
    gt_plasmid_length = ground_truth.groupby(by=0).max()[3]
    pred_plasmid_len = bin_file.groupby(by="plasmid").sum()["contig_len"]

    out += f"[GROUND TRUTH]\n\n"
    out += f"Total plasmids ground_truth length: {gt_plasmid_length} \n\n"

    out += f"[PREDICTION]\n\n"
    out += f"Total plasmids predicted length: {pred_plasmid_len} \n\n"

    ## calculate true positives, false negatives and false positives for SCORES
    t_pos = []
    f_neg = []
    f_pos = []
    t_pos_len = 0
    f_neg_len = 0
    f_pos_len = 0

    # gtruth_header = {0: label, 1: contig, 3:pl_lenght, 4:contig_len}
    for _, row in ground_truth.iterrows():
        if (
            row[0] != "chromosome"
        ):  ## here we can find TRUE POSITIVE and FALSE NEGATIVES plasmids
            if row[1] in bin_file["contig"].values:
                t_pos.append(row[1])
                t_pos_len += row[4]
            else:
                f_neg.append(row[1])
                f_neg_len += row[4]

        elif (
            row[0] == "chromosome"
        ):  ## here we can find FALSE PLASMIDS (FALSE POSITIVES)
            if row[1] in bin_file["contig"].values:
                f_pos.append(row[1])
                f_pos_len += row[4]

    out += f"True positives: {t_pos}\n"
    out += f"True positives LEN: {t_pos_len}\n"
    out += f"False negatives: {f_neg}\n"
    out += f"False negatives LEN: {f_neg_len}\n"
    out += f"False positives: {f_pos}\n"
    out += f"False positives LEN: {f_pos_len}\n\n"
    prec = safe_div(t_pos_len, (t_pos_len + f_pos_len))
    rec = safe_div(t_pos_len, (t_pos_len + f_neg_len))
    prec = safe_div(t_pos_len, (t_pos_len + f_pos_len))
    rec = safe_div(t_pos_len, (t_pos_len + f_neg_len))
    out += f"precision:\t\t{prec}\n"
    out += f"recall:\t\t\t{rec}\n"
    f1 = safe_div(2 * (prec * rec), (prec + rec))
    out += f"f1 score:\t\t{f1}\n"
    out += f"________________________________________\n"

    data = pd.DataFrame(
        {
            "sample": [sample] * 3,
            "type": ["f1_score", "precision", "recall"],
            "score": [f1, prec, rec],
        }
    )

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x="type",
            y=alt.Y("score", scale=alt.Scale(domain=(0, 1))),
            color=alt.Color("type").scale(range=["#e7ba52", "#c7c7c7", "#5f9ea0"]),
            order=alt.Order("type", sort="ascending"),
            column="sample",
        )
        .properties(title=f"[PRED_LABEL] (Plasmid SCORE")
    )
    chart.save(f"{args.output}.pred.scores.plasmids.pdf")

    if args.output:
        with open(f"{args.output}.stats.txt", "w", encoding="utf8") as f:
            f.write(out)

    else:
        print(out)


if __name__ == "__main__":
    main()
