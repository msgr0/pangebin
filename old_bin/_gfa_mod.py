"""
set of extension to use with gfapy. Try to include them in the gfapy package.

--remove-nodes
    input thr output
--rename
    input ren output

--pangenome_cleaner
    input output

        this file preprocess a pangenome gfa
        it adds the following tags to the segments:
        - OC: occurences in paths
        - cl: segment list
        - LN: length of the segment
        - aa: assembly type
        it adds the following tags to the paths:
        - LN: length of the path
        - aa: assembly type
            
        it also removes segments with length < 1 and segments with OC = 0
        it also removes paths that are not connected to any segment

--panassembly
    input assemblyA assemblyB output

"""
