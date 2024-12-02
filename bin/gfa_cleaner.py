import gfapy as gp
import argparse as ap


"""
Preprocess of a pangenome gfa to add the following tags to each segment:
- OC: occurences in paths
- cl: segment list
- LN: length of the segment
- aa: assembly type
it adds the following tags to the paths:
- LN: length of the path
- aa: assembly type
    
it also removes segments with length < 1 and segments with OC = 0
it also removes paths that are not connected to any segment
"""


def main():
    parser = ap.ArgumentParser(description="Clean a GFA file")
    parser.add_argument("--input", help="Input GFA file")
    parser.add_argument("--output", help="Output GFA file")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    gfa = gp.Gfa.from_file(input_file)
    gfa.validate()

    for seg in gfa.segments:
        seg.set_datatype("OC", "i")  # occurences in paths
        seg.OC = 0

        seg.set_datatype("cl", "Z")  # segment List
        seg.cl = ""

        seg.set_datatype("ll", "Z")
        seg.ll = ""

        seg.LN = len(seg.sequence)

        seg.set_datatype("aa", "A")  # assembler

    for path in gfa.paths:
        segs = path.segment_names
        _type = path.name[0]
        path_len = 0
        for seg in segs:
            segment = gfa.segment(seg.name)
            segment.OC += 1
            path_len += segment.LN
            if segment.cl == "":
                segment.cl = path.name
            else:
                segment.cl += "," + path.name

        path.set_datatype("LN", "i")
        path.LN = path_len

        path.set_datatype("aa", "A")
        path.aa = _type

    for path in gfa.paths:
        segs = path.segment_names
        for seg in segs:
            segment = gfa.segment(seg.name)
            segment_perc = segment.LN / float(path.LN)
            if segment.ll == "":
                segment.ll = f"{segment_perc:.6f}"
            else:
                segment.ll += f",{segment_perc:.6f}"




    for seg in gfa.segments:
        if seg.LN < 1:
            seg.disconnect()
            continue
        if seg.OC == 0:
            seg.disconnect()
            continue

        contig_list = seg.cl.split(",")
        contig_list = list(set(contig_list))
        sorted_clist = sorted(contig_list)
        if len(sorted_clist) >= 1:
            if sorted_clist[0] != "" and sorted_clist[0][0] != sorted_clist[-1][0]:
                seg.aa = "b"
            elif sorted_clist[0] != "":
                seg.aa = sorted_clist[0][0]

        

    gfa.validate()
    gfa.to_file(output_file)


if __name__ == "__main__":
    main()
