import sys
import argparse
import time

import entrez_utils as ez


def get_assembly_from_name(name: str):
    """Returns an assembly id from entrez"""
    assembly_id = ez.src("assembly", name)["IdList"]
    if len(assembly_id) == 0:
        return None
    elif len(assembly_id) >= 1:
        return assembly_id


def get_assembly_summary(assembly_id):
    """Returns an assembly summary object from entrez"""
    return ez.AssemblyObject(ez.summ("assembly", assembly_id))


def download_assembly(assembly_id):
    """Downloads an assembly from entrez"""
    try:
        summary = get_assembly_summary(assembly_id)
    except:
        summary = ""
        return None
    if len(summary.refseq_rsync) > 0:
        assembly_filename = summary.refseq_filename[:-3]
        tmp_rsync = summary.refseq_rsync
        tmp_file_gz = summary.refseq_filename
        return (tmp_rsync)
    elif len(summary.genebank_rsync) > 0:
        assembly_filename = summary.genebank_filename[:-3]
        tmp_rsync = summary.genebank_rsync
        tmp_file_gz = summary.genebank_filename
        return (tmp_rsync)
    else:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Accession ID to get the link assembly")

    args = parser.parse_args()
    assemblies = get_assembly_from_name(str(args.input))
    if assemblies is not None:
        ret_string = str(download_assembly(assemblies[0]))
        # print("------", file=sys.stderr)
        # print("sample:", str(args.input), file=sys.stderr)
        # print(str(download_assembly(assemblies[0])))
        # print("assembly:", assemblies[0], ret_string, file=sys.stderr)
        if ret_string == "None":
            pass
            # print("No assembly found", file=sys.stderr)
        else:
            print(ret_string, end="")
    else:
        pass
        # print("No assembly found", file=sys.stderr)


if __name__ == "__main__":
    main()
