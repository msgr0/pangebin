import sys
import argparse as ap

import pandas as pd
import gfapy
import os


"""
Enhance a pangenome graph built with PGGB. Requires the original Assembly graphs used to build the
pagenome graph in order to reconnect edges.
"""



### penalty_value = seg.LN/1000
### add Gfa_edges to Pangenome
### add Gfa_segments_tags to Pangenome Paths

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("--input",required=True, help="input pangenome graph")
    parser.add_argument("--assembly", required=True, help="gfa files of the assembly graphs given in input, with dp coverage annotated (unicycler formnat)", nargs="+")
    parser.add_argument("--output", required=True, help="output pangenome graph")
    
    args = parser.parse_args()
    assemblers = args.assembly
    pangenome = gfapy.Gfa.from_file(f"{args.input}")
    
    for a in assemblers:
        gfa = gfapy.Gfa.from_file(f"{a}")
        add_gfa_to_pangenome(gfa, pangenome)
    compute_metrics(pangenome)    
    pangenome.to_file(f"{args.output}")

    
##### GRAPH UTILS   

def get_path_by_name(gfa: gfapy.Gfa, name):
    """Get a path by name"""
    if gfa.line(str(name)) is not None:
        return gfa.line(str(name))
    else:
        print("path", name, "not found!!!!", type(name))
        return None
        # sys.exit(1)

def get_segment_by_name(gfa: gfapy.Gfa, name):
    """Get a path by name"""
    if gfa.segment(name) is not None:
        return gfa.segment(name)

def get_edge_by_def(gfa: gfapy.Gfa, def_: list):
    """Get a path by name"""
    for edge in gfa.dovetails:
        if def_[0] == edge.from_segment.name and def_[2] == edge.to_segment.name:
            if edge.from_orient == def_[1] and edge.to_orient == def_[3]:
                return edge


##### MAIN FUNCTIONS
        
def add_gfa_to_pangenome(gfa: gfapy.Gfa, pangenome: gfapy.Gfa):
    
    ### ADD GFA SEGMENT TAGS
    _type = gfa.segments[0].name[0]
    
    for seg in gfa.segments:
        
        ## segment in gfa is a path in the pangenome file
        match = pangenome.line(seg.name)
        if match is not None:
            assert match.aa == _type
            # if _type == "s":
            #     kmer_count = seg.KC
            #     assert seg.KC is not None
            #     match.set_datatype("KC", "i")
            #     match.KC = kmer_count
            if _type == "u" or _type == "s":
                dp = seg.dp
                assert seg.dp is not None
                match.set_datatype("dp", "f")
                match.dp = float(dp)
            else:
                raise TypeError("error! unknown assembler_type")
        
        else:
            print("segment", seg.name, "not found in pangenome")
            seg.disconnect()
            # assert False
            
    ### ADD GFA EDGES TO PANGENOME
    
    for edge in gfa.dovetails:
        
        from_contig = edge.from_segment
        to_contig = edge.to_segment

        from_orient = edge.from_orient
        to_orient = edge.to_orient
        
        try:
            path_from_contig = get_path_by_name(pangenome, from_contig.name)
        except:
            pass

        try:
            path_to_contig = get_path_by_name(pangenome, to_contig.name)
        except:
            pass

        if path_from_contig is None or path_to_contig is None:
            continue
        # link + + (last_fragment +) --- (first_fragment +)
        # link + - (last_fragment +) --- (last_fragment -)
        # link - + (first_fragment -) --- (first_fragment +)
        # link - - 1(first_fragment -) --- 2(last_fragment -) == 2( last_fragment +) --- 1( first_fragment +)

        if (from_orient == '+') and (to_orient == '+'):
            left_frag = [x.name for x in path_from_contig.segment_names][-1]
            left_orient = '+'
            right_frag = [x.name for x in path_to_contig.segment_names][0]
            right_orient = '+'

             
        elif (from_orient == '+') and (to_orient == '-'):
            left_frag = [x.name for x in path_from_contig.segment_names][-1]
            left_orient = '+'
            right_frag = [x.name for x in path_to_contig.segment_names][-1]
            right_orient = '-'


        elif (from_orient == '-') and (to_orient == '+'):
            left_frag = [x.name for x in path_from_contig.segment_names][0]
            left_orient = '-'
            right_frag = [x.name for x in path_to_contig.segment_names][0]
            right_orient = '+'

        elif (from_orient == '-') and (to_orient == '-'):
            #swapping left-right and using + + 
            left_frag = [x.name for x in path_to_contig.segment_names][-1]
            left_orient = '+'
            right_frag = [x.name for x in path_from_contig.segment_names][0]
            right_orient = '+'



        new_edge = gfapy.Line(
            f"L\t{left_frag}\t{left_orient}\t{right_frag}\t{right_orient}\t0M\taa:A:{_type}\tlt:Z:{_type}{_type}"
        )
    
        # end_fragment, end_fragment_orient = [(x.name, x.orient) for x in path_from_contig.segment_names][-1]
        # start_fragment, start_fragment_orient = [(x.name, x.orient) for x in path_to_contig.segment_names][0]
        
        
        # ## - strand
        # end_fragment_min, end_fragment_orient_min = [(x.name, x.orient) for x in path_from_contig.segment_names][0]
        # start_fragment_min, start_fragment_orient_min = [(x.name, x.orient) for x in path_to_contig.segment_names][-1]
        
        # if from_orient == "+":
        #     edge_begin = end_fragment
        #     edge_begin_orient = end_fragment_orient
        # elif from_orient == "-":
        #     edge_begin = end_fragment_min
        #     edge_begin_orient = end_fragment_orient_min
        # if to_orient == "+":
        #     edge_end = start_fragment
        #     edge_end_orient = start_fragment_orient
        # elif to_orient == "-":
        #     edge_end = start_fragment_min
        #     edge_end_orient = start_fragment_orient_min
            
        # new_edge = gfapy.Line(
        #     f"L\t{edge_begin}\t{edge_begin_orient}\t{edge_end}\t{edge_end_orient}\t0M\taa:A:{_type}\tlt:Z:{_type}{_type}"
        # )
        
        
        ##
        try:
            edge_to_update = get_edge_by_def(pangenome, [left_frag, left_orient, right_frag, right_orient])
            if edge_to_update is None:
                pangenome.add_line(new_edge)
            else:
                pass
                # print("edge", edge_to_update,"already present")
                # edge_to_update.set_datatype("ls", "A")
        except gfapy.error.NotUniqueError:
            pass
        

        
    
def compute_metrics(pangenome):
    
    for edge in pangenome.dovetails:
        if edge.lt is None:
            edge.set_datatype("aa", "A")
            edge.aa = "p"
            edge.set_datatype("lt", "Z")
            inc = pangenome.segment(edge.from_segment.name).aa
            out = pangenome.segment(edge.to_segment.name).aa
            if inc is None:
                inc = "n"
            if out is None:
                out = "n"
            edge.lt = inc + out
    
    #### annotate norm. coverage from Gfa Assemblies to Pangenome
    # total_coverage = 0
    # total_length = 0
    
    # for path in pangenome.paths:
    #     if path.aa == "s":
    #         total_coverage += path.KC
    #         total_length += path.LN
    
    for path in pangenome.paths:
        path.set_datatype("cv", "f")  # normalized Coverage
        try:
            assert path.dp is not None
            path.cv = path.dp
        except:
            print("IN PATH: path", path.name, "has no coverage:", path)
            sys.exit(1)


    # annotate mean coverage to fragments, based on paths (contigs) coverage
    for seg in pangenome.segments:
        seg.set_datatype("cv", "f")
        coverage_list = []
        if seg.cl is None:
            raise ValueError(f"segment {seg} has no paths associated")
        for contig in seg.cl.split(","):
            path = pangenome.line(contig)
            assert path is not None
            if contig == path.name:
                try:
                    assert path.cv is not None
                except:
                    print("IN SEG: path", path.name, "has no coverage:", path)
                    sys.exit(1)
                coverage_list.append(path.cv)
            else:
                raise KeyError(f"Path {contig} not found in pangenome")
        coverage_mean = sum([i for i in coverage_list])/len(coverage_list)
    
        assert coverage_mean is not None
        seg.cv = float(coverage_mean)

    ## annotate assembly penalty to segments
    for seg in pangenome.segments:
        seg.set_datatype("ap", "f")
        if False:
            segment.ap = 0
        else:
            if seg.aa == "u" or seg.aa == "s":
                penalty_value = seg.LN/1000
                penalty_value = 1 if penalty_value > 1 else penalty_value
                seg.ap = penalty_value
            else:
                seg.ap = 0

    


    
if __name__ == "__main__":
    main()
    sys.exit(0)


