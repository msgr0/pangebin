import gfapy as gp
import sys
import argparse as ap
from itertools import combinations
from itertools import product


true = True
false = False
null = None


"""
Remove GFA nodes based on a threshold length.
Reconnects edges of removed nodes with each respective Left Neighbour and Right Neighbour.
"""


def edge_exists(edge, collection):
    for e in collection:
        if edge_compare(edge, e) or edge_compare(edge_reverse(edge), e):
            return true
    return false


def edge_compare(edge1, edge2):
    if str(edge1[0][0]) == str(edge2[0][0]) and str(edge1[1][0]) == str(edge2[1][0]):
        if edge1[0][1] == edge2[0][1] and edge1[1][1] == edge2[1][1]:
            return true
    return false


def edge_reverse(edge):
    """
    from an edge( L(node1, orient1, position1), R(node2, orient2, position2))
    output the mirrored edge (R, L) with reversed orientation
    """
    assert edge[0][2] == "l", "edge[0][2] should be left, 'l'"
    assert edge[1][2] == "r", "edge[1][2] should be right, 'r'"

    return (
        (edge[1][0], invert(edge[1][1]), "l"),
        (edge[0][0], invert(edge[0][1]), "r"),
    )


def stoe(string):
    estr = str(string).split("\t")
    return ((estr[1], estr[2], "l"), (estr[3], estr[4], "r"))


def etos(edge):
    """
    from an edge( L(node1, orient1, position1), R(node2, orient2, position2))
    output the corresponding gfa formatted Link (L) string.
    """
    # edge[0] L_NODE
    # edge[0][0] name
    # edge[0][1] orient
    # edge[0][2] original pos (l=left, r=right)
    #
    # edge[1] R_NODE
    # edge[1][0] name
    # edge[1][1] orient
    # edge[1][2] original pos (l=left, r=right)
    #
    # if L_NODE(edge[0][2] == 'r') then reverse orient
    # if R_NODE(edge[1][2] == 'l') then reverse orient
    #
    # after that, if edge[0][1] == '-' == edge[1][1] then reverse polarity and invert position
    # i.e.    A -  B -   becomes   B +  A +    for clarity and consistency ...
    #
    l_orient = "?"
    r_orient = "?"

    if edge[0][2] == "r":
        l_orient = invert(edge[0][1])
    elif edge[0][2] == "l":
        l_orient = edge[0][1]

    if edge[1][2] == "l":
        r_orient = invert(edge[1][1])
    elif edge[1][2] == "r":
        r_orient = edge[1][1]

    return gp.Line(f"L\t{edge[0][0]}\t{l_orient}\t{edge[1][0]}\t{r_orient}\t0M")


def invert(sign):
    if sign == "+":
        return "-"
    elif sign == "-":
        return "+"
    if sign == "l":
        return "r"
    elif sign == "r":
        return "l"


def progress(counter, tot):
    if counter % 10 == 0:
        print(
            "working...",
            int(counter / tot * 10),
            "% graph analized",
            end="\r",
            file=sys.stderr,
        )
    return counter + 1


def extract_node(edge, orient, name):
    node = stoe(edge)

    if orient == "r":
        if node[0][0] == name and node[1][0] == name:
            return None
        elif node[1][1] == "-" and node[1][0] == name:
            # assert right[1][0] == seg.name
            # print("reversing R", node)
            node = edge_reverse(node)
            # print("into, ", node)
        elif node[0][1] == "+" and node[0][0] == name:
            assert node[0][0] == name
        else:
            assert False
        return (node[1][0], node[1][1], orient)

    elif orient == "l":
        if node[0][0] == name and node[1][0] == name:
            return None
        elif node[1][1] == "+" and node[1][0] == name:
            # assert right[1][0] == seg.name
            pass
        elif node[0][1] == "-" and node[0][0] == name:
            # print("reversing L", node)
            node = edge_reverse(node)
            # print("into, ", node)
        else:
            assert False
        return (node[0][0], node[0][1], orient)


def remove_lr(gfa, threshold):
    counter = 1
    for seg in gfa.segments:
        counter = progress(counter, len(gfa.segments))
        if (seg.LN != None and seg.LN <= threshold) or (len(seg.sequence) <= threshold):
            pass
        else:
            continue

        right_edges = list(seg.dovetails_R)
        left_edges = list(seg.dovetails_L)

        right_nodes = set()  # edge(node, orient, 'r')
        left_nodes = set()
        # print("removing seg", seg.name)
        # print("r-edge", right_edges)
        # print("l-edge", left_edges)

        for e in right_edges:
            n = extract_node(e, "r", seg.name)
            if n is not None:
                right_nodes.add(n)
                gfa.rm(e)
        for e in left_edges:
            n = extract_node(e, "l", seg.name)
            if n is not None:
                left_nodes.add(n)
                gfa.rm(e)

        # print("r-node", right_nodes)
        # print("l_node", left_nodes)

        gfa.rm(seg)
        # seg.disconnect()
        pairs = list(product(left_nodes, right_nodes))
        # print(pairs)
        gfa.validate()
        for edge in pairs:
            new_edge = etos(edge)
            # print (new_edge)
            try:
                gfa.add_line(new_edge)
                if edge[0][0] in ["13", "1"] and edge[1][0] in ["13", "1"]:
                    print("segment_removed:", seg.name, " providing edge:", new_edge)
            except:
                pass
                # print(f"edge already added!")
            # print(new_edge)
        gfa.validate()


def remove(gfa, threshold):

    counter = 0
    for seg in gfa.segments:
        counter = progress(counter, len(gfa.segments))

        if (seg.LN != None and seg.LN <= threshold) or (len(seg.sequence) <= threshold):
            nodes_to_reconnect = set()
            # print("segment: ", seg.name)
            for e in seg.dovetails:
                # print(e)
                # print(e.from_segment.name)
                if type(e.from_segment) == str:
                    continue
                if e.from_segment.name != seg.name:
                    nodes_to_reconnect.add(
                        (f"{e.from_segment.name}", f"{e.from_orient}", "l")
                    )
                elif e.to_segment.name != seg.name:
                    nodes_to_reconnect.add(
                        (f"{e.to_segment.name}", f"{e.to_orient}", "r")
                    )
                # else:
                #     try:
                #         gfa.rm(e)
                #     except:
                #         pass
                #     # print("Dup", e.from_segment.name, "--", e.to_segment.name)
                #     continue
                gfa.rm(e)  ## else do nothing and remove self edge.
                # print("---removed:", e)
            seg.disconnect()
            pairs = list(combinations(nodes_to_reconnect, 2))
            gfa.validate()
            if len(nodes_to_reconnect) == 1:
                # print("cappio")
                continue
            for edge in pairs:
                if edge[0][0] == edge[1][0]:
                    continue
                new_edge = etos(edge)
                # print (new_edge)
                try:
                    gfa.add_line(new_edge)
                    if edge[0][0] in ["13", "1"] and edge[1][0] in ["13", "1"]:
                        print(
                            "segment_removed:", seg.name, " providing edge:", new_edge
                        )

                except:
                    pass
                    # print(f"edge already added!")
                # print(new_edge)
            gfa.validate()


def main(args):
    gfa_file_path = args.input
    graph = gp.Gfa.from_file(gfa_file_path)

    if args.ver == "2":
        remove_lr(graph, args.threshold)
    elif args.ver == "1":
        remove(graph, args.threshold)
    else:
        print(
            "nothing done ... choose a version with --ver/-v flag = 1 for normal remove, = 2 for left-right removal"
        )
        sys.exit(1)
    graph.to_file(args.output)


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="remove gfa segments below the threshold and connect neighbours"
    )
    parser.add_argument("-i", "--input", help="input graph to trim")
    parser.add_argument("-o", "--output", help="output graph")
    parser.add_argument(
        "-t", "--threshold", help="threshold to remove contigs", type=int
    )
    parser.add_argument("-c", "--contains", help="string to remove contigs")
    parser.add_argument("-v", "--ver", help="version of the tool, either 1 or 2")

    args = parser.parse_args()

    main(args)
