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
import gfapy as gp
import os

def safe_div(x, y):
    """Safe division to handle division by zero"""
    if y == 0:
        return 0
    else:
        return (x / y)

def safe_div(x, y):
    """Safe division to handle division by zero"""
    if y == 0:
        return 0
    else:
        return (x / y)

def parser():
    
    parser = ap.ArgumentParser()
    parser.add_argument('--bin', help='Input bin file, tsv output of plasbin flow', required=True)
    parser.add_argument('--sample', help='Sample name', required=True)
    parser.add_argument('--csv', help='Input csv file, ground truth file', required=True)
    parser.add_argument('--gfa', help='Input gfa file containing the pangenome')
    parser.add_argument('--output', help='Output folder')
    parser.add_argument('--plot', help='Plot the results', action='store_true')
    parser.add_argument('--type', help='Print the results', required=True, choices=['pangenome', 'skesa', 'unicycler', 'assembler'])

    return parser.parse_args()

def plot_(data: pd.DataFrame, _output: str, _title: str):
        chart = alt.Chart(data).mark_bar().encode(
            x = 'sample',
            y = 'length',
            color = 'type',
            order = alt.Order('type', sort='ascending'),
            tooltip='length').properties(title=_title)
        chart.save(_output)
        chart.close()
        

def main():
    """ main function:
    - read the ground truth file
    - read the bin file, output of plasbin-flow
    - read the pangenome gfa file
    - calculate the statistics
    - print the results
    - plot the results
    """
    
    args = parser()
    
    sample = args.sample   
    out = f"Sample: {sample}\n"
    """GROUND TRUTH FILE [header]
    contig, label, length, plasmids, chromosomes, {assembler_type}"""
    ground_truth = pd.read_csv(args.csv)
    
    
    if args.type == "pangenome":
        """GROUND TRUTH TRANSFORMED DATAFRAME [header]"""
        header = ["contig","label","length","plasmids","chromosomes", "type"]
        gt_transformed = pd.DataFrame(columns=header)

        """PANGENOME [GFA]: transforming the ground truth file to the same format as the bin file
        i.e. in the pangenome fragments are part of the original contigs contained in the ground truth file"""
            
        gfa = gp.Gfa.from_file(args.gfa)
        gfa.validate()
        for seg in gfa.segments: # for each pangenome segment
            contig_list = seg.cl.split(",") ## retrieve the GT of each contig
            
            label = None
            
            for contig in contig_list:
                new_label = None

                if contig in ground_truth['contig'].values: # if the contig is in the ground truth file
                    new_label = ground_truth.loc[ground_truth['contig'] == contig, 'label'].values[0] # retrieve the label
                    if new_label in ["plasmid", "chromosome", "ambiguous"]: 
                        if label is None or label == "undef": # if this is the first label, assign it
                            label = new_label
                        elif label != new_label: # if it is not the first label, is ambiguous if it differs from the previous one
                            label = "ambiguous"
                        elif label == new_label: # (this assignment is not necessary, but it is more explicit)
                            label = new_label
                    else:
                        continue # do not annotate [short] and [undef] label
                else:
                    assert(False), f"Contig {contig} not found in the ground truth file"
                
                if new_label is None:
                    label = "undef"

            """Add the new segment(fragment) information to the transformed ground truth"""
            gt_transformed.loc[len(gt_transformed)] = {'contig': seg.name, 'label': label, 'length': seg.LN, 'type': seg.aa}
    else:
        """ OLD GT
        contig, label, length, plasmids, chromosomes, {assembler_type}"""

        """GROUND TRUTH TRANSFORMED DATAFRAME [header]"""
        header = ["contig","label","length","plasmids","chromosomes", "type"]
        gt_transformed = pd.DataFrame(columns=header)
        for _, row in ground_truth.iterrows():
            if row["contig"].startswith(str(args.type)[0:3]):
                if row["label"] in ["plasmid", "chromosome", "ambiguous"]:
                    gt_transformed.loc[len(gt_transformed)] = {'contig': row['contig'], 'label': row["label"], 'length': row['length'], 'type': str(args.type)[0]}
                else:
                    gt_transformed.loc[len(gt_transformed)] = {'contig': row['contig'], 'label': "undef", 'length': row['length'], 'type': str(args.type)[0]}
            
        

    """Compute the statistics, in len of fragments"""
    gt_plasmids_len = gt_transformed.groupby('label').sum()['length'].get('plasmid', 0)
    gt_chromosomes_len = gt_transformed.groupby('label').sum()['length'].get('chromosome',0)
    gt_ambigous_len = gt_transformed.groupby('label').sum()['length'].get('ambiguous', 0)
    gt_bd_plasmids = gt_plasmids_len + gt_ambigous_len
    gt_uncategorized = gt_transformed.groupby('label').sum()['length'].get('undef', 0)
    total_len = gt_bd_plasmids + gt_chromosomes_len + gt_uncategorized
    

    out += (f"\n")
    out += (f"[TYPE] {(args.type).upper()} statistics\n\n")
    addition = "(per assembler)" if args.type == "pangenome" else ""

    out += (f"[STATS] Total Contigs {addition}\n\n")
    if args.type == "pangenome":
        gt_contigs_len_u = gt_transformed.groupby('type').sum()['length'].get('u', 0 )
        gt_contigs_len_s = gt_transformed.groupby('type').sum()['length'].get('s', 0 )
        gt_contigs_len_b = gt_transformed.groupby('type').sum()['length'].get('b', 0)
        out += (f"\tUnicycler only (length): {gt_contigs_len_u}\n")
        out += (f"\tSkesa only (length): {gt_contigs_len_s}\n")
        out += (f"\tShared by both (length): {gt_contigs_len_b}\n")
        data = pd.DataFrame({"sample": [sample]*3, "type": ["both", "skesa_only", "unicyc_only"], "length": [gt_contigs_len_b, gt_contigs_len_s, gt_contigs_len_u]})
        plot_(data, os.path.join(args.output, f"{args.type}.stats.contigs.pdf"), f"[STAT] ({args.type})\n Contigs")
        
        
    out += (f"\t\t\t\t total: {total_len}\n\n")
    out += (f"________________________________________\n\n")

        

    out += (f"[LABEL] Chromosome {addition}\n\n")
    if args.type == "pangenome":
        gt_chrs_len_u = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('u', 'chromosome'), 0)
        # gt_chrs_len_u = gt_transformed.groupby(['type', 'label']).sum()['length'].get(['u', 'chromosome'], 0)
        gt_chrs_len_s = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('s', 'chromosome'), 0)
        gt_chrs_len_b = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('b', 'chromosome'), 0)
        out += (f"\tUnicycler only (length): {gt_chrs_len_u}\n")
        out += (f"\tSkesa only (length): {gt_chrs_len_s}\n")
        out += (f"\tShared by both (length): {gt_chrs_len_b}\n")
    out += (f"\t\t\t\t total: {gt_chromosomes_len}\n\n")
    out += (f"________________________________________\n\n")

    
    out += (f"[LABEL] Plasmid {addition}\n\n")
    if args.type == "pangenome" and "plasmid" in gt_transformed['label'].values:
            gt_plasmids_len_u = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('u', 'plasmid'), 0)
            gt_plasmids_len_s = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('s', 'plasmid'), 0)
            gt_plasmids_len_b = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('b', 'plasmid'), 0)
            out += (f"\tUnicycler only (length): {gt_plasmids_len_u}\n")
            out += (f"\tSkesa only (length): {gt_plasmids_len_s}\n")
            out += (f"\tShared by both (length): {gt_plasmids_len_b}\n")
            data = pd.DataFrame({"sample": [sample]*3, "type": ["both", "skesa_only", "unicyc_only"], "length": [gt_plasmids_len_b, gt_plasmids_len_s, gt_plasmids_len_u]})
            plot_(data, os.path.join(args.output, f"{args.type}.gt.plasmids.pdf"), f"[GT_LABEL] ({args.type})\n Plasmid")
            
    out += (f"\t\t\t\t total: {gt_plasmids_len}\n\n")

    out += (f"________________________________________\n\n")

    out += (f"[LABEL] Plasmid + Ambiguous {addition}\n\n")
    if args.type == "pangenome" and "ambiguous" in gt_transformed['label'].values:
        gt_ambigous_len_u = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('u', 'ambiguous'), 0)
        gt_ambigous_len_s = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('s', 'ambiguous'), 0)
        gt_ambigous_len_b = gt_transformed.groupby(['type', 'label']).sum()['length'].get(('b', 'ambiguous'), 0)
        
        gt_bd_plasmids_len_u = gt_plasmids_len_u + gt_ambigous_len_u
        gt_bd_plasmids_len_s = gt_plasmids_len_s + gt_ambigous_len_s
        gt_bd_plasmids_len_b = gt_plasmids_len_b + gt_ambigous_len_b
        
        out += (f"\tUnicycler only (length): {gt_bd_plasmids_len_u}\n")
        out += (f"\tSkesa only (length): {gt_bd_plasmids_len_s}\n")
        out += (f"\tShared by both (length): {gt_bd_plasmids_len_b}\n")
        
        data = pd.DataFrame({"sample": [sample]*3, "type": ["both", "skesa_only", "unicyc_only"], "length": [gt_ambigous_len_b, gt_ambigous_len_s, gt_ambigous_len_u]})
        plot_(data, os.path.join(args.output, f"{args.type}.gt.ambiguous.pdf"), f"[GT_LABEL] ({args.type})\n Ambiguous")
                
        data = pd.DataFrame({"sample": [sample]*3, "type": ["both", "skesa_only", "unicyc_only"], "length": [gt_bd_plasmids_len_b, gt_bd_plasmids_len_s, gt_bd_plasmids_len_u]})
        plot_(data, os.path.join(args.output, f"{args.type}.gt.broad_plasmids.pdf"), f"[GT_LABEL] ({args.type})\n Plasmid + Ambiguous")
        
    out += (f"\t\t\t\t total: {gt_bd_plasmids} \n\n")
    out += (f"________________________________________\n\n")



    # Read in the bin file
    """#Pls_ID	Flow	GC_bin	Length	Contigs"""
    bin_file = pd.read_csv(args.bin, sep='\t')

    header = ["contig", "label", "length","plasmids", "type"]
    header = ["contig", "label", "length","plasmids", "type"]
    bin_transformed = pd.DataFrame(columns=header)
    ## if bin_file has no rows
    if bin_file.empty:
        out += (f"[PREDICTION]\n\n")
        out += (f"Total plasmids predicted length: 0 \n\n")
        out += (f"________________________________________\n")        
        if args.output:
            with open(f"{args.output}/{args.type}.stats.txt", "w", encoding='utf8') as f:
                f.write(out)
        else:
            print(out)
        return
    
    for _, row in bin_file.iterrows():
        contigs = [x.split(":")[0] for x in row['Contigs'].split(',')]
        for contig in contigs:
            type_ = gfa.segment(contig).aa if args.type == "pangenome" else str(args.type)[0]
            bin_transformed.loc[len(bin_transformed)] = {'contig': contig, 'label': "plasmid", 'plasmids': row['#Pls_ID'], 'length': gt_transformed.loc[gt_transformed['contig'] == contig, 'length'].values[0], 'type': type_}
            bin_transformed.loc[len(bin_transformed)] = {'contig': contig, 'label': "plasmid", 'plasmids': row['#Pls_ID'], 'length': gt_transformed.loc[gt_transformed['contig'] == contig, 'length'].values[0], 'type': type_}
    
    bin_transformed.sort_values(by='contig', inplace=True)
    # remove duplicates
    bin_transformed.drop_duplicates(subset='contig', keep='first', inplace=True)
    pred_plasmids_len = bin_transformed.groupby('label').sum()['length']['plasmid']
    bin_transformed.sort_values(by='contig', inplace=True)
    # remove duplicates
    bin_transformed.drop_duplicates(subset='contig', keep='first', inplace=True)
    pred_plasmids_len = bin_transformed.groupby('label').sum()['length']['plasmid']
    other_len = total_len - pred_plasmids_len
    
    
    
    out += (f"[PREDICTION]\n\n")
    out += (f"Total plasmids predicted length: {pred_plasmids_len} \n\n")
    
    ## calculate true positives, false negatives and false positives for SCORES
    t_pos = {'b':[], 'u':[], 's':[]}
    f_neg = {'b':[], 'u':[], 's':[]}
    f_pos = {'b':[], 'u':[], 's':[]}
    t_pos_len = {'b':0, 'u':0, 's':0}
    f_neg_len = {'b':0, 'u':0, 's':0}
    f_pos_len = {'b':0, 'u':0, 's':0}
    
    for _, row in gt_transformed.iterrows():
        if row['label'] == 'plasmid' or row['label'] == 'ambiguous': ## here we can find TRUE POSITIVE and FALSE NEGATIVES plasmids
            if row['contig'] in bin_transformed['contig'].values:
                t_pos[row['type']].append(row['contig'])
                t_pos_len[row['type']] += row['length']
            else:
                f_neg[row['type']].append(row['contig'])
                f_neg_len[row['type']] += row['length']

        elif row['label'] == 'chromosome': ## here we can find FALSE PLASMIDS (FALSE POSITIVES)
            if row['contig'] in bin_transformed['contig'].values:
                f_pos[row['type']].append(row['contig'])
                f_pos_len[row['type']] += row['length']

    if args.type == "pangenome":
        out+=(f"Skesa\n")
        out+=(f":---True positives: ({t_pos_len['s']}) {t_pos['s']}\n")
        out+=(f":---False negatives: ({f_neg_len['s']}) {f_neg['s']}\n")
        out+=(f":---False positives: ({f_pos_len['s']}) {f_pos['s']}\n\n")
        skesa_precision = safe_div(t_pos_len['s'], (t_pos_len['s'] + f_pos_len['s']))
        skesa_recall = safe_div(t_pos_len['s'], (t_pos_len['s'] + f_neg_len['s']))
        out+=(f"precision:\t\t{skesa_precision}\n")
        out+=(f"recall:\t\t\t{skesa_recall}\n\n")
        
        
        out+=(f"Unicycler\n")
        out+=(f":---True positives: ({t_pos_len['u']}) {t_pos['u']}\n")
        out+=(f":---False negatives: ({f_neg_len['u']}) {f_neg['u']}\n")
        out+=(f":---False positives: ({f_pos_len['u']}) {f_pos['u']}\n\n")
        uni_precision = safe_div(t_pos_len['u'], (t_pos_len['u'] + f_pos_len['u']))
        uni_recall = safe_div(t_pos_len['u'], (t_pos_len['u'] + f_neg_len['u']))
        out+=(f"precision:\t\t{uni_precision}\n")
        out+=(f"recall:\t\t\t{uni_recall}\n\n")
        
        out+=(f"Both\n")
        out+=(f":---True positives: ({t_pos_len['b']}) {t_pos['b']}\n")
        out+=(f":---False negatives: ({f_neg_len['b']}) {f_neg['b']}\n")
        out+=(f":---False positives: ({f_pos_len['b']}) {f_pos['b']}\n\n")
        both_precision = safe_div(t_pos_len['b'], (t_pos_len['b'] + f_pos_len['b']))
        both_recall = safe_div(t_pos_len['b'], (t_pos_len['b'] + f_neg_len['b']))
        out+=(f"precision:\t\t{both_precision}\n")
        out+=(f"recall:\t\t\t{both_recall}\n")
        
        data = pd.DataFrame({"sample": [sample]*3, "type": ["both", "skesa_only", "unicyc_only"], "length": [t_pos_len['b']+f_pos_len['b'],
                                                                                                             t_pos_len['s']+f_pos_len['s'],
                                                                                                             t_pos_len['u']+f_pos_len['u']]})
        
        plot_(data, os.path.join(args.output, f"{args.type}.pred.plasmids.pdf"), f"[PRED_LABEL] ({args.type})\n Plasmids")                                                                                                     
        

        
        total_precision = safe_div((t_pos_len['s'] + t_pos_len['u'] + t_pos_len['b']), (t_pos_len['s'] + t_pos_len['u'] + t_pos_len['b'] + f_pos_len['s'] + f_pos_len['u'] + f_pos_len['b']))
        total_recall = safe_div((t_pos_len['s'] + t_pos_len['u'] + t_pos_len['b']), (t_pos_len['s'] + t_pos_len['u'] + t_pos_len['b'] + f_neg_len['s'] + f_neg_len['u'] + f_neg_len['b']))
        total_f1_score = safe_div(2 * (total_precision * total_recall), (total_precision + total_recall))
        
        ### precision plot
        data = pd.DataFrame({"sample": [sample]*4, "type": ["both", "skesa_only", "unicyc_only", "z_final"], "score": [both_precision,
                                                                                                             skesa_precision,
                                                                                                             uni_precision,
                                                                                                             total_precision]})
        chart = alt.Chart(data).mark_bar().encode(
            x='type',
            y=alt.Y('score', scale=alt.Scale(domain=(0, 1))),
            color=alt.Color('type'),
            order=alt.Order(
                'type',
                sort='ascending'
            ),
            column='sample'
            
        ).properties(
            title=f"[PRED_LABEL] ({args.type})\n Plasmid PRECISION"
        )
        chart.save(os.path.join(args.output, f"{args.type}.pred.precision.plasmids.pdf"))
        
        # recall plot
        
        
        data = pd.DataFrame({"sample": [sample]*4, "type": ["both", "skesa_only", "unicyc_only", "z_final"], "score": [both_recall,
                                                                                                                skesa_recall,
                                                                                                                uni_recall,
                                                                                                                total_recall]})
        
        chart = alt.Chart(data).mark_bar().encode(
            x='type',
            y=alt.Y('score', scale=alt.Scale(domain=(0, 1))),
            color=alt.Color('type'),
            order=alt.Order(
                'type',
                sort='ascending'
            ),
            column='sample'
            
        ).properties(
            title=f"[PRED_LABEL] ({args.type})\n Plasmid RECALL"
        )
        chart.save(os.path.join(args.output, f"{args.type}.pred.recall.plasmids.pdf"))  
                
        
        ## scores plot
        data = pd.DataFrame({"sample": [sample]*3, "type": ["f1_score", "precision", "recall"], "score": [total_f1_score, total_precision, total_recall]})
                
        chart = alt.Chart(data).mark_bar().encode(
            x='type',
            y=alt.Y('score', scale=alt.Scale(domain=(0, 1))),
            color=alt.Color('type').scale(range=['#e7ba52', '#c7c7c7', '#5f9ea0']),
            order=alt.Order(
                'type',
                sort='ascending'
            ),
            column='sample'
            
        ).properties(
            title=f"[PRED_LABEL] ({args.type})\n Plasmid SCORE"
        )
        chart.save(os.path.join(args.output, f"{args.type}.pred.scores.plasmids.pdf"))
        

        out += (f"________________________________________\n")        
        out += (f"Overall Statistics\n\n")
        out += (f"precision:\t\t{total_precision}\n")
        out += (f"recall:\t\t\t{total_recall}\n")
        out += (f"F1 Score:\t\t{total_f1_score}\n")
        out += (f"________________________________________\n")        
        
   
        
    else:
        acc = str(args.type[0])
        out+=(f"{str(args.type).capitalize()}\n")

        out+=(f"True positives: {t_pos[acc]}\n")
        out+=(f"False negatives: {f_neg[acc]}\n")
        out+=(f"False positives: {f_pos[acc]}\n\n")
        prec = safe_div(t_pos_len[acc],(t_pos_len[acc] + f_pos_len[acc]))
        rec = safe_div(t_pos_len[acc], (t_pos_len[acc] + f_neg_len[acc]))
        prec = safe_div(t_pos_len[acc],(t_pos_len[acc] + f_pos_len[acc]))
        rec = safe_div(t_pos_len[acc], (t_pos_len[acc] + f_neg_len[acc]))
        out+=(f"precision:\t\t{prec}\n")
        out+=(f"recall:\t\t\t{rec}\n")
        f1 = safe_div(2 * (prec * rec), (prec + rec))
        out+=(f"f1 score:\t\t{f1}\n")
        out += (f"________________________________________\n")
        
        data = pd.DataFrame({"sample": [sample]*3, "type": ["f1_score", "precision", "recall"], "score": [f1, prec, rec]})
        
        chart = alt.Chart(data).mark_bar().encode(
            x='type',
            y=alt.Y('score', scale=alt.Scale(domain=(0, 1))),
            color=alt.Color('type').scale(range=['#e7ba52', '#c7c7c7', '#5f9ea0']),
            order=alt.Order(
                'type',
                sort='ascending'
            ),
            column='sample'
            
        ).properties(
            title=f"[PRED_LABEL] ({args.type})\n Plasmid SCORE"
        )
        chart.save(os.path.join(args.output, f"{args.type}.pred.scores.plasmids.pdf"))
        


            
    
    if args.output:
        with open(f"{args.output}/{args.type}.stats.txt", "w", encoding='utf8') as f:
            f.write(out)
            
    else:
        print(out)
    
    gt_transformed.to_csv(f"{args.output}/{args.type}.transformed.gt.csv", sep=',', index=False, encoding='utf-8')
    bin_transformed.to_csv(f"{args.output}/{args.type}.transformed.bins.csv", sep=',', index=False, encoding='utf-8')
    mix_transformed_header = ["contig", "label", "pred_label", "length", "type", "plasmids", "chromosomes"]
    
    mix_transformed = pd.DataFrame(columns=mix_transformed_header)
    mix_transformed = pd.concat([mix_transformed, gt_transformed], ignore_index=True)

    for _, row in bin_transformed.iterrows():
        if row['contig'] in mix_transformed['contig'].values:
            mix_transformed.loc[mix_transformed['contig'] == row['contig'], 'pred_label'] = row["label"]       
    mix_transformed.to_csv(f"{args.output}/{args.type}.transformed.mix.csv", sep=',', index=False, encoding='utf-8')
    gt_transformed.to_csv(f"{args.output}/{args.type}.transformed.gt.csv", sep=',', index=False, encoding='utf-8')
    bin_transformed.to_csv(f"{args.output}/{args.type}.transformed.bins.csv", sep=',', index=False, encoding='utf-8')
    mix_transformed_header = ["contig", "label", "pred_label", "length", "type", "plasmids", "chromosomes"]
    
    mix_transformed = pd.DataFrame(columns=mix_transformed_header)
    mix_transformed = pd.concat([mix_transformed, gt_transformed], ignore_index=True)

    for _, row in bin_transformed.iterrows():
        if row['contig'] in mix_transformed['contig'].values:
            mix_transformed.loc[mix_transformed['contig'] == row['contig'], 'pred_label'] = row["label"]       
    mix_transformed.to_csv(f"{args.output}/{args.type}.transformed.mix.csv", sep=',', index=False, encoding='utf-8')
if __name__ == "__main__":
    main()
    
    
