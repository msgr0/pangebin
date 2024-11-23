import gfapy as gp

import argparse as ap

def main():
    parser = ap.ArgumentParser(description="Rename a GFA file with _type")
    parser.add_argument("--input", help="Input GFA file")
    parser.add_argument("--output", help="Output GFA file")
    parser.add_argument("--type", help="Type of the GFA file")
    
    args = parser.parse_args()  
    
    graph = gp.Gfa.from_file(args.input, vlevel=0)
            
    out = rename_gfa(graph, args.type)
    out.vlevel = 3
    out.validate()
    out.to_file(args.output)
    

def rename_gfa(gfa, _type):
    gfa.vlevel = 3
    counter = 1
    for seg in gfa.segments:
        if seg.LN is not None:
            if seg.LN == 0:
                seg.sequence = "*"
                seg.LN = 0
                gfa.validate()
        else:
            if len(seg.sequence) == 0:
                seg.sequence = "*"
                seg.LN = 0
                gfa.validate()
        seg.name = f"{_type}_{counter}"
        counter += 1
    # for e in gfa.edges:
    #     if gfa.segment(e.from_segment) is None:
    #         print(e.from_segment)
    #     elif gfa.segment(e.to_segment) is None:
    #         print(e.to_segment)

    return gfa
        
if __name__ == "__main__":
    main()
    