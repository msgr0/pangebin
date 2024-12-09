# This module is used to eval plasmid bins, output of plasbinflow and gplasi
# tools.

# JUST LABELING MODE!!!

# MODE COMPARE
# input: PREDICTION1, PREDICTION2, GROUNDTRUTH(previously built using
# build_truth_*.py tools)
# --graph out_graph with the scores, comparing PREDICTION1 and PREDICTION2

# MODE EVALUATE
# input: PREDICTION, GROUNDTRUTH(previously built using build_truth_*.py tools)
# --graph out_graph with the scores of PREDICTION against the GT

import argparse as ap

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Evaluate one prediction files OR Compare two prediction files against a Ground Truth File")
    # select either "compare" or "evaluate"
    parser.add_argument("--mode", "-m", "mode selected for the tool")


def evaluate(pred, gt):
    graph()
    pass


def compare(pred_a, pred_b, gt):
    pass


def graph():
    pass
