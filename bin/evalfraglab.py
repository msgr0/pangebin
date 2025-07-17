"""
EVALuate FRAGments LABeling (evalfraglab)
========================================

This script evaluates the performance of pangebin on the pan-assembly graph
by comparing predicted plasmid contigs against plasmid contigs computed on skesa and unicycler assemblies.
It reads predictions and ground truth data from TSV files, calculates precision, recall, and F1
score, and outputs the results to a specified file.


Usage:
    evalfraglab.py <predictions_file> <ground_truth_file> <output_file>

Example:
    python evalfraglab.py  ground_truth.tsv evaluation_results.txt
"""

import argparse
import pandas as pd
import sys, os, logging
import gfapy as gp
from typing import Dict, List, Tuple

def parse_args():
    parser = argparse.ArgumentParser(description="Compare pangebin binning vs unicycler and skesa binning.")
    parser.add_argument('-d', '--dir', help="Path to dataset dir.")
    return parser.parse_args()


def main():
    
    args = parse_args()
    data = args.dir
    header_gt = [
        "plasmid_label",
        "contig",
        "plasmid_length",
        "contig_length",
        "coverage",
    ]

    header_bin = [
        "plasmid",
        "contig",
        "contig_len",
    ]
    print("starting") 
    # Read the pangebin predictions
    

    sample_perf = pd.DataFrame(columns=["sample", "metric", "TP", "FP", "FN", "prec", "rec", "f1"])
    perf = {"TP":0, "FP":0, "FN":0, "precision":0.0, "recall":0.0, "f1":0.0}
    unip = perf.copy()
    skep = perf.copy()
    panp = perf.copy()
    
    for sample in os.listdir(data):
        # sample_perf[sample] = perf.copy()
        
        if sample.startswith('.'):
            continue
        if not os.path.isdir(os.path.join(data, sample)):
            continue
        print("----------Processing file:", sample)
        fske = os.path.join(data, sample, f"{sample}.ske.1.pbf.pred.tab")
        funi = os.path.join(data, sample, f"{sample}.uni.1.pbf.pred.tab")
        fpan = os.path.join(data, sample, f"{sample}.pan.1.pbf.pred.tab")
        fsgt = os.path.join(data, sample, f"{sample}.1.ske.gt.tsv") 
        fugt = os.path.join(data, sample, f"{sample}.1.uni.gt.tsv") 
        fmgt = os.path.join(data, sample, f"{sample}.1.mix.pan.gt.tsv")
        fpsg = os.path.join(data, sample, f"{sample}.1.ske.pan.gt.tsv")
        fpug = os.path.join(data, sample, f"{sample}.1.uni.pan.gt.tsv")

        panasm = os.path.join(data, sample, f"{sample}.panasm.gfa")
        cf, fc, allc, allf = mapgfa(panasm)
        
        skedf = pd.read_csv(fske, sep="\t", header=0, dtype={'plasmid':str, 'contig':str, 'contig_len':str})
        unidf = pd.read_csv(funi, sep="\t", header=0, dtype={'plasmid':str, 'contig':str, 'contig_len':str})
        pandf = pd.read_csv(fpan, sep="\t", header=0, dtype={'plasmid':str, 'contig':str, 'contig_len':str})

        gt_cols = {'plasmids':str, 'contig':str, 'score':float, 'total_len':int, 'contig_len':int}

        gtudf = pd.read_csv(fugt, sep="\t", header=None)
        gtudf.columns = gt_cols
        gtsdf = pd.read_csv(fsgt, sep="\t", header=None)
        gtsdf.columns = gt_cols
        gtpmdf = pd.read_csv(fmgt, sep="\t", header=None)
        gtpmdf.columns = gt_cols
        gtpudf = pd.read_csv(fpug, sep="\t", header=None)
        gtpudf.columns = gt_cols
        gtpsdf = pd.read_csv(fpsg, sep="\t", header=None)
        gtpsdf.columns = gt_cols
        #convert contig column to string
        gtudf['contig'] = gtudf['contig'].astype(str)
        gtsdf['contig'] = gtsdf['contig'].astype(str)
        gtpmdf['contig'] = gtpmdf['contig'].astype(str)
        gtpudf['contig'] = gtpudf['contig'].astype(str)
        gtpsdf['contig'] = gtpsdf['contig'].astype(str)

        if pandf.empty:
            print(f"Skipping {sample} due to empty pan predictions")
            continue 
        
        skedf['fragment'] = skedf['contig'].apply(lambda x: cf.get(x,x))
        unidf['fragment'] = unidf['contig'].apply(lambda x: cf.get(x,x))
        pandf['fragment'] = pandf['contig'].apply(lambda x: fc.get(x,x))
        

        #remove chromosome from gt
        
        gtudf['fragment'] = gtudf['contig'].apply(lambda x: cf.get(x,x))
        gtsdf['fragment'] = gtsdf['contig'].apply(lambda x: cf.get(x,x))
        gtpmdf['fragment'] = gtpmdf['contig'].apply(lambda x: fc.get(x,x))
        gtpudf['fragment'] = gtpudf['contig'].apply(lambda x: fc.get(x,x))
        gtpsdf['fragment'] = gtpsdf['contig'].apply(lambda x: fc.get(x,x))

        len_dict = {}
        
        for _, row in gtudf.iterrows():
            contig = row['contig']
            if contig not in len_dict:
                len_dict[contig] = row['contig_len']
        for _, row in gtsdf.iterrows():
            contig = row['contig']
            if contig not in len_dict:
                len_dict[contig] = row['contig_len']
        for _, row in gtpmdf.iterrows():
            contig = row['contig']
            if contig not in len_dict:
                len_dict[contig] = row['contig_len']
                
        gtudf = gtudf[~gtudf['plasmids'].str.startswith('chromosome')]
        gtsdf = gtsdf[~gtsdf['plasmids'].str.startswith('chromosome')]
        gtpmdf = gtpmdf[~gtpmdf['plasmids'].str.startswith('chromosome')]
        gtpudf = gtpudf[~gtpudf['plasmids'].str.startswith('chromosome')]
        gtpsdf = gtpsdf[~gtpsdf['plasmids'].str.startswith('chromosome')]

        pred_pan = dict()
        pred_pan['tp'], pred_pan['fp'], pred_pan['fn'] = eval_pred(pandf, gtpmdf, len_dict)
        pred_pan = comp_metrics(pred_pan['tp'], pred_pan['fp'], pred_pan['fn'])
        pred_panS = dict()
        pred_panS['tp'], pred_panS['fp'], pred_panS['fn'] = eval_pred(pandf, gtpsdf, len_dict)
        pred_panS = comp_metrics(pred_panS['tp'], pred_panS['fp'], pred_panS['fn'])
        pred_panU = dict()
        pred_panU['tp'], pred_panU['fp'], pred_panU['fn'] = eval_pred(pandf, gtpudf, len_dict)
        pred_panU = comp_metrics(pred_panU['tp'], pred_panU['fp'], pred_panU['fn'])
        pred_ske = dict()
        pred_ske['tp'], pred_ske['fp'], pred_ske['fn'] = eval_pred(skedf, gtsdf, len_dict)
        pred_ske = comp_metrics(pred_ske['tp'], pred_ske['fp'], pred_ske['fn'])
        pred_uni = dict()
        pred_uni['tp'], pred_uni['fp'], pred_uni['fn'] = eval_pred(unidf, gtudf, len_dict)
        pred_uni = comp_metrics(pred_uni['tp'], pred_uni['fp'], pred_uni['fn'])
        
        
        print("pangenome mix\t", pred_pan)
        print("pangenome skesa\t", pred_panS)
        print("pangenome uni\t", pred_panU)
        print("skesa\t\t", pred_ske)
        print("unicycler\t", pred_uni)
        
        pan_ske_diff = diff_pred(pandf, skedf, gtpsdf, len_dict)
        print("pan_ske diff\t", pan_ske_diff)

        pan_uni_diff = diff_pred(pandf, unidf, gtpudf, len_dict)
        print("pan_uni diff\t", pan_uni_diff)

        pan_ass_diff = diff_pred(pandf, pd.concat([skedf, unidf]), gtpmdf, len_dict)
        pan_ass_perc = eval_perc(pandf, gtpmdf, len_dict)
        print("pan_ass diff\t", pan_ass_diff)
        print(pan_ass_diff['TP%'] if pan_ass_diff['TP%'] > 1.0 else "not")
        
        # add pangenome_mix as a row of sample_perf
        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
            "sample": sample,
            "metric": "pangenome_mix",
            "TP": pred_pan["TP"],
            "FP": pred_pan["FP"],
            "FN": pred_pan["FN"],
            "prec": pred_pan["prec"],
            "rec": pred_pan["rec"],
            "f1": pred_pan["f1"]
            }])
        ], ignore_index=True)

        # add every other set
        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
            "sample": sample,
            "metric": "pangenome_skes",
            "TP": pred_panS["TP"],
            "FP": pred_panS["FP"],
            "FN": pred_panS["FN"],
            "prec": pred_panS["prec"],
            "rec": pred_panS["rec"],
            "f1": pred_panS["f1"]
            }])
        ], ignore_index=True)

        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
            "sample": sample,
            "metric": "pangenome_uni",
            "TP": pred_panU["TP"],
            "FP": pred_panU["FP"],
            "FN": pred_panU["FN"],
            "prec": pred_panU["prec"],
            "rec": pred_panU["rec"],
            "f1": pred_panU["f1"]
            }])
        ], ignore_index=True)

        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
            "sample": sample,
            "metric": "skesa",
            "TP": pred_ske["TP"],
            "FP": pred_ske["FP"],
            "FN": pred_ske["FN"],
            "prec": pred_ske["prec"],
            "rec": pred_ske["rec"],
            "f1": pred_ske["f1"]
            }])
        ], ignore_index=True)

        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
            "sample": sample,
            "metric": "unicycler",
            "TP": pred_uni["TP"],
            "FP": pred_uni["FP"],
            "FN": pred_uni["FN"],
            "prec": pred_uni["prec"],
            "rec": pred_uni["rec"],
            "f1": pred_uni["f1"]
            }])
        ], ignore_index=True)

        # add pan_ske_diff as a row of sample_perf
        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
                "sample": sample,
                "metric": "pan_ske_diff",
                "TP": pan_ske_diff['TP+'],
                "FP": pan_ske_diff['FP+'],
                "FN": pan_ske_diff['FN+'],
                "prec": pan_ske_diff['prec+'],
                "rec": pan_ske_diff['rec+'],
                "f1": pan_ske_diff['f1+']
            }])
        ], ignore_index=True)
        # add pan_uni_diff as a row of sample_perf
        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
                "sample": sample,
                "metric": "pan_uni_diff",
                "TP": pan_uni_diff['TP+'],
                "FP": pan_uni_diff['FP+'],
                "FN": pan_uni_diff['FN+'],
                "prec": pan_uni_diff['prec+'],
                "rec": pan_uni_diff['rec+'],
                "f1": pan_uni_diff['f1+']
            }])
        ], ignore_index=True)
        # #to csv sample perf
        sample_perf = pd.concat([
            sample_perf,
            pd.DataFrame([{
                "sample": sample,
                "metric": "pan_ass_diff",
                "TP": pan_ass_diff['TP+'],
                "FP": pan_ass_diff['FP+'],
                "FN": pan_ass_diff['FN+'],
                "prec": pan_ass_diff['prec+'],
                "rec": pan_ass_diff['rec+'],
                "f1": pan_ass_diff['f1+']
            }])
        ], ignore_index=True)


    # compute mean on sample_perf, per metric (across samples)
    mean_perf = sample_perf.groupby("metric")[["TP", "FP", "FN", "prec", "rec", "f1"]].mean().reset_index()
    sample_perf.to_csv(f"sample_perf.csv", index=False)
    mean_perf.to_csv(f"mean_perf.csv", index=False)
    



def diff_pred(pred_df: pd.DataFrame, diff_df: pd.DataFrame, gt_df: pd.DataFrame, lendict: dict) -> pd.DataFrame:
    """
    Computes the difference between predicted and ground truth fragments.
    Returns a DataFrame with the differences.
    """
    pred_frags = set(pred_df['contig'])
    #expand fragment column of diff_frags into strings
    # diff_frags = set(diff_df['fragments'])
    #extract fragments from diff_df
    diff_frags = set([x for sublist in diff_df['fragment'] for x in sublist])
    rem_frags = pred_frags - diff_frags
    TP = cumlen(rem_frags & set(gt_df['contig']), lendict)
    FP = cumlen(rem_frags - set(gt_df['contig']), lendict)
    FN = cumlen(set(gt_df['contig']) - rem_frags, lendict)
    TP_PERC = 100 * TP / cumlen(set(gt_df['contig']), lendict) if set(gt_df['contig']) else 0.0
    FP_PERC = 100 * FP / cumlen(pred_frags, lendict) if pred_frags else 0.0
    FN_PERC = 100 * FN / cumlen(set(gt_df['contig']), lendict) if set(gt_df['contig']) else 0.0
    
    prec = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    rec = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * (prec * rec) / (prec + rec) if (prec + rec) > 0 else 0.0
    
    


    return {'TP%': float(f"{TP_PERC:.3f}"), 'FP%': float(f"{FP_PERC:.3f}"), 'FN%': float(f"{FN_PERC:.3f}"), 'TP+': TP, 'FP+': FP, 'FN+': FN, 'prec+': float(f"{prec:.3f}"), 'rec+': float(f"{rec:.3f}"), 'f1+': float(f"{f1:.3f}")}

    # diff_frags = pred_frags - diff_df['contig']
    # diff_df = pd.DataFrame(list(diff_frags), columns=['contig'])
    
    # # Add lengths from lendict
    # diff_df['length'] = diff_df['contig'].apply(lambda x: lendict.get(x, 0))
    
    # return diff_df
        

def eval_pred(pred_df: pd.DataFrame, gt_df: pd.DataFrame, lendict: dict) -> Tuple[int, int, int]:  
    pred_frags = set(pred_df['contig'])
    gt_frags = set(gt_df['contig'])
    tp = cumlen(pred_frags & gt_frags, lendict)
    fp = cumlen(pred_frags - gt_frags, lendict)
    fn = cumlen(gt_frags - pred_frags, lendict)
    return tp, fp, fn

def eval_perc(pred_df: pd.DataFrame, gt_df: pd.DataFrame, lendict: dict) -> float:
    pred_frags = set(pred_df['contig'])
    gt_frags = set(gt_df['contig'])

    tp = cumlen(pred_frags & gt_frags, lendict)
    return tp / cumlen(gt_frags, lendict) if gt_frags else 0.0


def comp_metrics(tp: int, fp: int, fn: int) -> Dict[str, float]:
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    # 
    return {"TP": tp, "FP": fp, "FN": fn, "prec": float(f"{precision:.3f}"), "rec": float(f"{recall:.3f}"), "f1": float(f"{f1:.3f}")}

def cumlen(frags: set, lenmap: dict) -> int:
    return sum(lenmap.get(frag,0) for frag in frags)

def mapgfa(gfaf: str):
    cf = dict()
    fc = dict()
    all_ctg = set()
    all_frag = set()
    
    graph = gp.Gfa.from_file(gfaf)
    if not graph:
       raise ValueError(f"Failed to read GFA file: {gfaf}")
    
    for segment in graph.segments:
        fragment_id = segment.name
        all_frag.add(fragment_id)
        fc[fragment_id] = set()
        fragments_list = set(str(segment.cl).split(','))
        for contig in fragments_list:
            all_ctg.add(contig)
            cf.setdefault(contig, set()).add(fragment_id)
        fc[fragment_id].update(fragments_list)

    return cf, fc, all_ctg, all_frag

if __name__ == "__main__":
    main()