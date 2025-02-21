"""Gets the fna link from entrez, given an assembly entry from STDIN"""

import sys
import argparse as ap
import time
import subprocess

import entrez_utils as ez

def main(args):
  accession = args.input
  
  assemblies = get_assembly_from_name(str(accession))
  if assemblies is not None:
    ret_string = str(download_assembly(assemblies[0]))
    if ret_string == "None":
      # print('0', end='')
      pass
    else:
      subprocess.call(["curl", ret_string, "--output", args.output])
      # rsync -vPp --chmod=u+w \$link ${genomic_fna_gz}
      # print("")
      # print(ret_string, end="")
  else:
    pass
    # print('-1', end='')
  

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

if __name__ == "__main__":
  parser = ap.ArgumentParser(description="ncbi genome retrival script")
  parser.add_argument("--input")
  parser.add_argument("--output")
  args = parser.parse_args()
  
  main(args)
