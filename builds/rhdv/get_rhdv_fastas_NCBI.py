from Bio import Entrez
from Bio import SeqIO

# download RHDV genomes from NCBI with required header metadata
# Matthew J. Neave 28.06.2018

# get NCBI sequence IDs from Mahar et al. 2017. J. Virol
# these ids range from MF421563 - MF421701

mahar_ids = ["MF" + str(rec) for rec in range(421563, 421702)]

# use biopython go online and retrieve these records
# need to enter an email so users don't overuse their servers

Entrez.email = "matthewjneave1@gmail.com"

def retrieve_NCBI_record(NCBI_ID):
    new_handle = Entrez.efetch(db="nucleotide", id=NCBI_ID, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    collection_date = seq_record.features[0].qualifiers["collection_date"]
    print(seq_record.features[0].qualifiers)



retrieve_NCBI_record("MF421603")