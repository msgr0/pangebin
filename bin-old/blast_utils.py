import argparse
import pandas as pd
import logging as log
import sys

def main():
    log.basicConfig(filename='blast_rename.log', encoding='utf-8', level=log.DEBUG)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", "--assembly", help="Path to assembly fna file", required=True
    )
    parser.add_argument(
        "-c", "--csv", help="Path to contigs list", required=True
    )
    parser.add_argument(
        "-f", "--fasta", help="Path to sequences list", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Path to output file", required=False
    )
        
    args = parser.parse_args()
    
    plasmid_count = 0
    chromosome_count = 0
    
    contigs = []
    contigs_len = []
    contigs_tags = []
    curr_len = 0 
    with open(args.fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                contigs_len.append(curr_len)
                curr_len = 0
                contigs.append(line.split()[0][1:])
                contigs_tags.append(line.split()[1:])
            else:
                curr_len += len(str(line))
    contigs_len.append(curr_len)                
    contigs_len = contigs_len[1:]
    
    log.debug(f"Found {len(contigs)} contigs")
    log.debug(f"Example contig: {contigs[0]}, length: {contigs_len[0]}, tags: {contigs_tags[0]}")  
             
    assembly = []
    assembly_desc = []
    assembly_len = []
    curr_len = 0
    with open(args.assembly, "r") as f:
        for line in f:
            if line.startswith(">"):
                assembly_len.append(curr_len)                
                curr_len = 0
                assembly.append(line.split()[0][1:])
                assembly_desc.append(" ".join(line.split()[1:]))
            else:
                curr_len += len(str(line))  
    assembly_len.append(curr_len)                
    assembly_len = assembly_len[1:]
    
    log.debug(f"Found {len(assembly)} sequences")

    blast_sample = pd.read_csv(
            args.csv, sep="\t", header=None
        )
    
    for i in range(0, len(assembly)):
        description = [x.replace(",", "") for x in assembly_desc[i].lower().split()]

        id = assembly[i]
        
        if any(ele in description for ele in ["chromosome", "chromosomes"]):
            log.debug(f"Found chromosome: {id}")
            chromosome_count += 1
            blast_sample.loc[blast_sample[1] == id, 12] = "chromosome"
            blast_sample.loc[blast_sample[1] == id, 13] = f"chr_{chromosome_count}"    
        elif any(ele in description for ele in ["plasmid", "plasmids"]):
            log.debug(f"Found plasmid: {id}")
            plasmid_count += 1
            blast_sample.loc[blast_sample[1] == id, 12] = "plasmid"
            blast_sample.loc[blast_sample[1] == id, 13] = f"pl_{plasmid_count}"
        
    for i in range(0, len(contigs)):
        id = contigs[i]
        length = contigs_len[i]
        blast_sample.loc[blast_sample[0] == id, 14] = length

    
    if plasmid_count == 0:
        log.debug(f"No plasmids found")
        print(f"No plasmids found")
        sys.exit(1)
        
        
    out_cols = [
        "contig",
        "label",
        "length",
        "plasmids",
        "chromosomes",
    ]
    mixed = pd.DataFrame(columns=out_cols)
    
    """
    blast_sample header
    
    qseqid  sseqid          pident  length  mismatch qstart qend    bitscore
    ske_1	NZ_CP024859.1	100.000	20058	0	      1	    20058	37041
    """
    blast_sample = blast_sample.groupby(0)
    for name, group in blast_sample:
        save_line = True
        chromosomes = []
        plasmids = []
        for index, row in group.iterrows():
            if row[3] <= 50 :
                group.drop(index=index, inplace=True)
                continue
            if row[2] < 95.0:
                group.drop(index=index, inplace=True)
                continue

        if len(group.nunique(0)) > 0:
            # those are ambiguous contigs
            for index, row in group.iterrows():
                contig_type = row[12]
                if row[12] == "chromosome":
                    if row[3] >= 0.90 * row[14]:
                        chromosomes.append(row[13])
                    elif row[7] >= 0.95 * row[14]:
                        chromosomes.append(row[13])
                elif row[12] == "plasmid":
                    if row[3] >= 0.90 * row[14]:
                        plasmids.append(row[13])
                    elif row[7] >= 0.95 * row[14]:
                        plasmids.append(row[13])
            if len(chromosomes) > 0 and len(plasmids) > 0:
                contig_type = "ambiguous"
            elif len(chromosomes) > 0:
                contig_type = "chromosome"
            elif len(plasmids) > 0:
                contig_type = "plasmid"
            else:
                contig_type = "short"

            # if contig_type != "short":
            
            new_row = {
                "contig": name,
                "label": contig_type,
                "length": group[14].max(),
                "plasmids": ','.join(set(plasmids)),
                "chromosomes": ','.join(set(chromosomes)),
            }

            if contig_type != "short":   
                mixed.loc[len(mixed)] = new_row

    for c in contigs:
        if str(c) not in mixed['contig'].values:
            label = "short" if contigs_len[contigs.index(c)] < 50 else "undef"
            new_row = {
                "contig": c,
                "label": label,
                "length": contigs_len[contigs.index(c)],
                "plasmids": "",
                "chromosomes": "",
            }
            mixed.loc[len(mixed)] = new_row
    pd.DataFrame(mixed).to_csv(
        f"{args.output}",
        header=mixed.columns,
        index=None,
        sep=",",
        na_rep="",
        mode="w",
    )


    pd.DataFrame(blast_sample).to_csv(
        f"blast_map.tsv",
        header=None,
        index=None,
        sep="\t",
        mode="w",
    )

        
if __name__ == "__main__":
    main()  

