from Bio import Entrez
from Bio import SeqIO

# download RHDV genomes from NCBI with required header metadata
# Matthew J. Neave 02.07.2018

# get NCBI sequence IDs from the oldest RHDV2 in Aus

# use biopython go online and retrieve these records
# need to enter an email so users don't overload their servers

Entrez.email = "matthewjneave1@gmail.com"


def retrieve_ncbi_record(ncbi_id):
    new_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    # extracting several bits of information from different parts of the record
    collection_date = "2015-05-13"
    isolate = seq_record.features[0].qualifiers["isolate"][0]
    state = "ACT"
    country = "Australia"
    genotype = seq_record.features[0].qualifiers["note"][0].lstrip("genotype: ").replace(" ", "_")
    authors = seq_record.annotations['references'][0].authors.split(",")[0] + " et al"
    title = seq_record.annotations['references'][0].title
    journal = seq_record.annotations['references'][0].journal
    header = ">" + "|".join([isolate, "RHDV", seq_record.id, genotype, collection_date, country, state,
                             authors,
                             title, journal])
    return(header, seq_record)


# use this function to retrieve all the records and write to fasta file
# with headers correctly formatted for Nextstrain

oldest_fasta = open("oldest_RHDV.fasta", "w")

record = retrieve_ncbi_record("KT280060")
oldest_fasta.write(record[0] + "\n" + str(record[1].seq) + "\n")


