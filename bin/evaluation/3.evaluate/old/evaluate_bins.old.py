import argparse as ap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import altair as alt
import sklearn.metrics as skm


"""
take a bin file (tsv) output of plasbin flow and the"""


def main():
    parser = ap.ArgumentParser()
    parser.add_argument("--bin", help="Input bin file, tsv output of plasbin flow")
    parser.add_argument("--csv", help="Input csv file, ground truth file")
    parser.add_argument("--output", help="Output folder")
    parser.add_argument("--plot", help="Plot the results", action="store_true")
    parser.add_argument(
        "--type", help='Type of bins, either (p)"pangenome" or (a)"assembler"'
    )

    args = parser.parse_args()

    ## Read the ground truth file
    "contig,label,length,plasmids,chromosomes"
    ground_truth = pd.read_csv(args.csv)

    gt_plasmids_len = ground_truth.groupby("label").sum()["length"]["plasmid"]
    gt_chromosomes_len = ground_truth.groupby("label").sum()["length"]["chromosome"]
    gt_undef_len = ground_truth.groupby("label").sum()["length"]["undef"]
    gt_short_len = ground_truth.groupby("label").sum()["length"]["short"]
    gt_ambigous_len = ground_truth.groupby("label").sum()["length"]["ambiguous"]

    gt_plasmid = gt_plasmids_len + gt_ambigous_len
    gt_other = gt_chromosomes_len + gt_undef_len + gt_short_len

    # Read in the bin file
    "#Pls_ID	Flow	GC_bin	Length	Contigs"
    bin_file = pd.read_csv(args.bin, sep="\t")

    header = ["contig", "pred", "length", "plasmids"]
    bin_transformed = pd.DataFrame(columns=header)
    for index, row in bin_file.iterrows():
        contigs = [x.split(":")[0] for x in row["Contigs"].split(",")]
        for contig in contigs:
            bin_transformed.loc[len(bin_transformed)] = {
                "contig": contig,
                "pred": "plasmid",
                "plasmids": row["#Pls_ID"],
                "length": ground_truth.loc[
                    ground_truth["contig"] == contig, "length"
                ].values[0],
            }

    # print('Predicted')
    # print('---- count ----')
    # print(bin_transformed.groupby('pred').count()['contig'])
    # print('---- cumulated length ----')
    # print(bin_transformed.groupby('pred').sum()['length'])

    pred_plasmids_len = bin_transformed.groupby("pred").sum()["length"]["plasmid"]
    other_len = gt_other + gt_plasmid - pred_plasmids_len

    OUT_ = f"######## SAMPLE ########\n"
    OUT_ += f"\n"
    OUT_ += f"#### GT Statistics ########\n"
    OUT_ += f"Total assembly length: {gt_plasmid + gt_other} \n"
    OUT_ += f"######## GT statistics ##########\n"
    OUT_ += f"Of which plasmids (length): {gt_plasmids_len}\n"
    OUT_ += f"Of which ambigous (length): {gt_ambigous_len}\n"
    OUT_ += f"Total plasmid + ambiguous (length): {gt_plasmid}\n"
    OUT_ += f"________________________________________\n\n"

    OUT_ += f"######## Prediction Statistics ########\n"
    OUT_ += f"Total assembly prediction length (double_check): {pred_plasmids_len + other_len} \n"
    OUT_ += f"Of which plasmids (length): {pred_plasmids_len}\n"
    OUT_ += f"________________________________________\n\n"

    ## calculate true positives, false negatives and false positives for SCORES

    true_positives = []
    false_negatives = []
    false_positives = []

    true_positives_len = 0
    false_negatives_len = 0
    false_positives_len = 0

    for index, row in ground_truth.iterrows():
        if (
            row["label"] == "plasmid" or row["label"] == "ambiguous"
        ):  ## here we can find TRUE POSITIVE and FALSE NEGATIVES plasmids
            if row["contig"] in bin_transformed["contig"].values:
                true_positives.append(row["contig"])
                true_positives_len += row["length"]
            else:
                false_negatives.append(row["contig"]).get("cavolo", "ciao")
                false_negatives_len += row["length"].get("hello", pippo)

        elif (
            row["label"] == "chromosome"
        ):  ## here we can find FALSE PLASMIDS (FALSE POSITIVES)
            if row["contig"] in bin_transformed["contig"].values:
                false_positives.append(row["contig"])
                false_positives_len += row["length"]

    # We are discarding the short and undef contigs, as they are not relevant for the evaluation

    OUT_ += f"######## Prediction statistics ##########\n"
    OUT_ += f"True positives length: {true_positives_len}\n"
    OUT_ += f"False negatives length: {false_negatives_len}\n"
    OUT_ += f"False positives length: {false_positives_len}\n"
    OUT_ += f"________________________________________\n"

    try:
        precision = true_positives_len / (true_positives_len + false_positives_len)
    except ZeroDivisionError:
        precision = 0

    try:
        recall = true_positives_len / (true_positives_len + false_negatives_len)
    except ZeroDivisionError:
        recall = 0

    try:
        f1_score = 2 * (precision * recall) / (precision + recall)
    except ZeroDivisionError:
        f1_score = 0

    OUT_ += f"######## CUMULATIVE Prediction statistics ##########\n"
    OUT_ += f"PRECISION\n"
    OUT_ += f"{precision} \n"
    OUT_ += f"________________________________________\n"
    OUT_ += f"RECALL\n"
    OUT_ += f"{recall}\n"
    OUT_ += f"________________________________________\n"
    OUT_ += f"F1 Score\n"
    OUT_ += f"{f1_score}\n"
    OUT_ += f"############################################\n\n"

    print(OUT_)


if __name__ == "__main__":
    main()
