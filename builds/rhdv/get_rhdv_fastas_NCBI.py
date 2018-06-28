from Bio import Entrez
from Bio import SeqIO

# download RHDV genomes from NCBI with required header metadata
# Matthew J. Neave 28.06.2018

# get NCBI sequence IDs from Mahar et al. 2017. J. Virol
# these ids range from MF421563 - MF421701

mahar_ids = ["MF" + str(rec) for rec in range(421563, 421702)]

# use biopython go online and retrieve these records
# need to enter an email so users don't overload their servers

Entrez.email = "matthewjneave1@gmail.com"

def retrieve_NCBI_record(NCBI_ID):
    new_handle = Entrez.efetch(db="nucleotide", id=NCBI_ID, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    # extracting several bits of information from different parts of the record
    collection_date = seq_record.features[0].qualifiers["collection_date"][0]
    isolate = seq_record.features[0].qualifiers["isolate"][0]
    state = isolate.split("/")[1]
    country = seq_record.features[0].qualifiers["country"][0]
    genotype = seq_record.features[0].qualifiers["note"][0].lstrip("genotype: ").replace(" ", "_")
    authors = seq_record.annotations['references'][0].authors.split(",")[0] + " et al"
    title = seq_record.annotations['references'][0].title
    journal = seq_record.annotations['references'][0].journal
    citation = " ".join([title, journal])
    header = ">" + "|".join([isolate, "RHDV", seq_record.id, genotype, collection_date, country, state, authors, citation])
    return(header, seq_record)

# use this function to retrieve all the records and write to fasta file

mahar_fasta = open("mahar_RHDV.fasta", "w")

for record_id in mahar_ids:
    record = retrieve_NCBI_record(record_id)
    mahar_fasta.write(record[0] + "\n" + str(record[1].seq) + "\n")


