from Bio import Entrez
from Bio import SeqIO

# download RHDV genomes from NCBI with required header metadata

record_id = "MF421563"

Entrez.email = "matthewjneave1@gmail.com"

new_handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="genbank")

seq_record = SeqIO.read(new_handle, "genbank")

print(seq_record.features[0].qualifiers["collection_date"])