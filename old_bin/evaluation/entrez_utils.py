""" Module to hold functions for interacting with entrez """
from Bio import Entrez
import time

# BEGIN parameters
Entrez.email = "mttsgro@gmail.com"
Entrez.tool = "plasbin-flow@pangenome"
Entrez.api_key = "7c5ce3b0327c240c68dc660082564a0bda08"
WAIT_TIME = 0.5


class AssemblyObject:
    """Class to hold information about an assembly retrived from entrez"""

    def __init__(self, summary):
        root = summary["DocumentSummarySet"]["DocumentSummary"][0]
        self.organism = root["Organism"]
        self.tax_id = root["Taxid"]
        self.refseq_id = root["RsUid"]
        self.genebank_id = root["GbUid"]
        self.biosample = root["BioSampleAccn"]
        self.genebank_acc = root["Synonym"]["Genbank"]
        self.refseq_acc = root["Synonym"]["RefSeq"]
        self.assembly_acc = root["AssemblyAccession"]
        self.assembly_name = root["AssemblyName"]
        self.sample_name = root["BioSampleAccn"]
        self.ftp_genebank = root["FtpPath_GenBank"]
        self.ftp_refseq = root["FtpPath_RefSeq"]
        self.genebank_genomic = (
            self.ftp_genebank
            + "/"
            + self.genebank_acc
            + "_"
            + self.assembly_name
            + "_genomic.fna.gz"
            if len(self.ftp_genebank) > 0
            else ""
        )
        self.refseq_genomic = (
            self.ftp_refseq
            + "/"
            + self.refseq_acc
            + "_"
            + self.assembly_name
            + "_genomic.fna.gz"
            if len(self.ftp_refseq) > 0
            else ""
        )
        self.genebank_rsync = (
            "https://" + self.genebank_genomic.split("//")[1]
            if len(self.ftp_genebank) > 0
            else ""
        )
        self.refseq_rsync = (
            "https://" + self.refseq_genomic.split("//")[1]
            if len(self.ftp_refseq) > 0
            else ""
        )
        self.genebank_filename = (
            self.genebank_genomic.split("/")[-1] if len(self.ftp_genebank) > 0 else ""
        )
        self.refseq_filename = (
            self.refseq_genomic.split("/")[-1] if len(self.ftp_refseq) > 0 else ""
        )

    def __str__(self):
        return self.assembly_name


def src(database, term, retmax=10000):
    """Searches entrez for a term and returns the results"""

    time.sleep(WAIT_TIME)
    handle = Entrez.esearch(db=database, retmax=retmax, term=term)
    record = Entrez.read(handle, validate=False)
    handle.close()
    return record


def summ(database, id, retmax=10000):
    """Retrieves a summary of an assembly from entrez"""
    time.sleep(WAIT_TIME)
    handle = Entrez.esummary(db=database, retmax=retmax, id=id)
    record = Entrez.read(handle, validate=False)
    handle.close()
    return record
