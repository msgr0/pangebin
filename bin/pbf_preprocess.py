
import gfapy as gp
import pandas as pd
import argparse as ap
import sys



"""
Normalize gd and gc scores from plasbinflow preprocess by the conti list in the pagenome
"""
def main():
    parser = ap.ArgumentParser()
    parser.add_argument("--input",required=True, help="input pangenome graph")
    parser.add_argument("--gc", required=False, help="tsv-scores-files for GC content from pbf-preprocess of the assembly graphs given in input", nargs="+")
    parser.add_argument("--gd", required=False, help="tsv-scores-files for GDensity from pbf-preprocess of the assembly graphs given in input", nargs="+")
    parser.add_argument("--output", required=True, help="output  pangenome graph")
    
    args = parser.parse_args()

    pangenome = gp.Gfa.from_file(f"{args.input}")

    if (not args.input or not args.output):
        print("enter input output parameters")
        sys.exit(0)
    if (args.gd): # modify plasmidness
        # gd_file pangenome_fragment score ...
        scores = pd.read_csv(f"{args.gd}", sep='\t', header=None, usecols=[0,1], index_col=0, names=['contig', 'score'])
        p_scores = pd.DataFrame(columns=['contig', 'score'], index=0)
        for seg in pangenome.segments:
            score_sum = 0.0
            for contig, perc_len in seg.cl.split(','), seg.ll.split(','):
                score_sum += perc_len * scores.loc[contig, 'score']
                perc_sum += perc_len
            score_sum /= perc_len
            p_scores[seg] = score_sum
        p_scores.to_csv(f"{a}")
    # if (args.gc): # modify gc-content
    #     # gd_file pangenome_fragment score score ...
    #     for seg in pangenome.segments:
    #         scores_sum = []
    #         for contig, perc_len in seg.cl.split(','), seg.ll.split(','):
    #             perc_len = 0
    #             for i in range(len(scores_sum)):
    #                 scores_sum[i] += perc_len * score[i](contig)
    #                 perc_sum += perc_len
    #         for i in range(score_sum):
    #             score_sum[i] /= perc_len
    #         tsv_out[seg] = score_sum[i]

    if (not args.gc and not args.gd):
        print("give at least one among --gd and --gd tsv files")
        sys.exit(0)




