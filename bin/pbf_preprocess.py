
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
        if (not args.input and not args.gd):
            print("enter at least one file to convert, gene_density(gd) or gc_content(gc)")

        sys.exit(0)

    if (args.gd): # modify plasmidness
        # gd_file pangenome_fragment score ...
        for gd_file in args.gd:
            print("... converting gene density file(s) ...")
            scores = pd.read_csv(f"{gd_file}", index_col=0, header=None,  usecols=[0,1], sep='\t')
            # scores.to_csv("test_in.csv", sep='\t', header=None)
            p_scores = pd.DataFrame(columns=['score'])
            # p_scores.to_csv('test_out.csv', sep='\t', header=None)

            for seg in pangenome.segments:
                score_sum = 0.0
                len_perc_sum = 0.0
        
                contigs = str(seg.cl).split(',')
                len_percs = str(seg.ll).split(',')
                if len(len_percs) == 0 or len(contigs) == 0:
                    continue
                # print(contigs, len_percs)
                
                for l in len_percs:
                    len_perc_sum += float(l)
                
                for c, l in zip(contigs, len_percs):
                    score_sum += float(l) * scores.loc[c, 1]

                score_sum /= float(len_perc_sum) if len_perc_sum > 0.0 else 1.0
                # print(seg.name, score_sum)

                p_scores.loc[seg.name] = score_sum
    
            p_scores.to_csv(f"{args.output}.gd.tsv", sep='\t', header=None)
    if (args.gc): # modify gc-content
        for gc_file in args.gc:
            print("... converting gc content file(s) ...")
            gc = pd.read_csv(f"{gc_file}", index_col=0, header=None, sep='\t')
            gc_scores = pd.DataFrame(columns=[1,2,3,4,5,6])

            for seg in pangenome.segments:
                scores = [0 for x in list(gc.columns)]
                len_perc_sum = 0.0
                contigs = str(seg.cl).split(',')
                len_percs = str(seg.ll).split(',')
                if len(len_percs) == 0 or len(contigs) == 0:
                    continue
                # print(contigs, len_percs)
                
                for l in len_percs:
                    len_perc_sum += float(l)
                # print(len_perc_sum)

                for c, l in zip(contigs, len_percs):
                    for x in range(len(scores)):
                        scores[x] += float(l) * gc.at[c, x+1]
                        # print(scores[x])

                for x in range(len(scores)):
                    scores[x] = (scores[x] / float(len_perc_sum)) if len_perc_sum > 0.0 else scores[x]
                    gc_scores.at[seg.name, x+1] = scores[x]
                    # print(seg.name, scores[x])

            gc_scores.to_csv(f"{args.output}.gc.tsv", sep='\t', header=None)

    if (not args.gc and not args.gd):
        print("give at least one among --gd and --gd tsv files")
        sys.exit(0)



if __name__ == "__main__":
    main()
