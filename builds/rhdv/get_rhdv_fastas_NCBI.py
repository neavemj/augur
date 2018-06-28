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
    collection_date = seq_record.features[0].qualifiers["collection_date"][0]
    isolate = seq_record.features[0].qualifiers["isolate"][0]
    state = isolate.split("/")[1]
    country = seq_record.features[0].qualifiers["country"][0]
    genotype = seq_record.features[0].qualifiers["note"][0].lstrip("genotype: ").replace(" ", "_")
    print(seq_record.annotations['references'][0])
    authors = seq_record.annotations['references'][0].authors.split(",")[0] + " et al"
    title = seq_record.annotations['references'][0].title
    journal = seq_record.annotations['references'][0].journal
    citation = title +
    header = ">" + "|".join([isolate, "RHDV", seq_record.id, genotype, collection_date, country, state, authors])
    print(header)


retrieve_NCBI_record("MF421603")