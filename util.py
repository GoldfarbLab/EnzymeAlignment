def fasta_id_to_filename(id):
    return id.replace("|", "_")

def fasta_id_to_fasta_filename(id):
    return fasta_id_to_filename(id) + ".fasta"
