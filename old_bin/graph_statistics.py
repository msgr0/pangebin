import argparse
import gfapy as gp

"""
Provide graphs statistics
"""


def main(args):
    gfa_path = args.input
    if args.pangenome:
        print("outputting statiscs in *pangenome* mode")
    else:
        print("outputting statics in graph mode")

    gfa = gp.Gfa.from_file(gfa_path)

    node_count = 0
    edges_count_1 = 0
    edges_count_2 = 0
    node_degree = {}
    thresholds = [10, 25, 50, 100, 250, 500, 1000, 5000]
    node_lengths = {x: 0 for x in thresholds}
    node_lengths.update({"over": 0})
    node_cumulative_l = {x: 0 for x in thresholds}
    node_cumulative_l.update({"over": 0})

    for seg in gfa.segments:
        node_count += 1
        seg_deg = len(set([x.name for x in seg.neighbours]))
        edges_count_1 += seg_deg

        try:
            node_degree[seg_deg] += 1
        except:
            node_degree.update({seg_deg: 1})

        node_len = seg.LN if seg.LN != None else len(seg.sequence)
        for x in thresholds:
            if node_len < x:
                node_lengths[x] += 1
                break
            elif x == thresholds[-1]:
                node_lengths["over"] += 1

    for x in node_cumulative_l:
        for y in thresholds:
            if x == "over":
                node_cumulative_l[x] += node_lengths[y]
                continue
            elif x >= y:
                node_cumulative_l[x] += node_lengths[y]
        if x == "over":
            node_cumulative_l[x] += node_lengths[x]

    edges_count_1 *= 0.5  # counting the neighbours of each node shopuld result in double the edge count.
    if args.pangenome:
        edges_type = {"p": 0, "u": 0, "s": 0}
    for line in gfa.lines:
        if str(line).startswith("L"):
            edges_count_2 += 1
            if line.aa != None:
                edges_type[line.aa] += 1

    print("node_count\n", node_count)
    print("edges_count_1\n", edges_count_1)
    print("edges_count_2\n", edges_count_2)

    node_degree_s = {i: node_degree[i] for i in sorted(node_degree.keys())}
    print("node_degree\n", node_degree_s)

    # print("node_degree\n", sort(node_degree.keys()))
    print("thresholds\n", thresholds)
    print("node_lengths\n", node_lengths)
    print("node_cumulative_lengths{len: count}\n", node_cumulative_l)
    if args.pangenome:
        print("edges per type", edges_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute simple graph statistcs")
    parser.add_argument("-i", "--input", required=True, help="input file .gfa")
    parser.add_argument(
        "-p",
        "--pangenome",
        action="store_true",
        help="if is a pan-assembly graph, output individual assembly statistics",
    )

    args = parser.parse_args()
    main(args)
