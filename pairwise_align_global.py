#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO, AlignIO
from math import floor, ceil
from seed_protein import SeedProtein
import util
import time
import tempfile
import shutil


########################################################################################################################
# Functions
########################################################################################################################


def align_all_seeds(seed_id2proteins, candidate, out_path, job_id, min_alignment_identity_threshold):

    candidate_prefix_path = os.path.join(out_path, str(job_id), util.fasta_id_to_filename(candidate.id))
    candidate_fasta_path = os.path.join(out_path, str(job_id), util.fasta_id_to_fasta_filename(candidate.id))

    candidate_fasta_written = False

    for seed in seed_id2proteins.values():
        if (len(candidate.seq) * min_alignment_identity_threshold > len(seed.seq) or  # too long
                len(candidate.seq) < len(seed.seq) * min_alignment_identity_threshold):  # too short
            continue

        # if the candidate passes the sequence length thresholds, only write the candidate fasta file once
        if not candidate_fasta_written:
            #SeqIO.write(candidate, candidate_fasta_path, "fasta")
            candidate_fasta_written = True

        align_to_seed(seed, ">" + candidate.id + "\n" + str(candidate.seq), candidate_prefix_path, candidate_fasta_path, min_alignment_identity_threshold)


def align_to_seed(seed, candidate_str, candidate_prefix_path, candidate_fasta_path, min_alignment_identity_threshold):

    aln_path = candidate_prefix_path + "_vs_" + util.fasta_id_to_filename(seed.id) + ".aln.txt"

    with tempfile.SpooledTemporaryFile(max_size=1e6, mode="w+") as aln_file:
        p = subprocess.Popen(["stretcher",
                         "-bsequence", seed.fasta_path,
                         # "-asequence", seed.fasta_path,
                         # "-bsequence", candidate_fasta_path,
                         # "-outfile", aln_path,
                         "-gapopen", "1",
                         "-gapextend", "1",
                         "-aformat3", "pair",
                         "-auto",
                         "-filter",
                         "-sprotein1",
                         "-sprotein2"], stdout=aln_file, stdin=subprocess.PIPE, universal_newlines=True)

        p.communicate(input=candidate_str)

        #print(candidate_fasta_path)

        # check if the candidate passes the sequence identity threshold
        aln_file.seek(0)
        for line in aln_file:
            if line.startswith("# Identity: "):
                identity = float(line.rstrip().split("(")[-1][0:-2]) / 100.0
                #print(identity)
                if identity < min_alignment_identity_threshold:
                    return
                break

        # check for conserved residues
        aln_file.seek(0)
        alignment = AlignIO.read(aln_file, "emboss")

        seed_seq_index = 0
        candidate_seq_index = 0
        last_conserved_index_match = 0

        for i in range(0, len(alignment[0].seq)):
            if alignment[1].seq[i] != "-":
                candidate_seq_index += 1
            if alignment[0].seq[i] != "-":
                seed_seq_index += 1
                if seed_seq_index in seed.index2conserved_residue:
                    conserved_residue = seed.index2conserved_residue[seed_seq_index]

                    # Find closest match in allowable subsequence. Tie breaks go to the left.
                    min_index = max(0, i - conserved_residue.n_term_tolerance, last_conserved_index_match)
                    max_index = min(len(alignment[0].seq)-1, i + conserved_residue.c_term_tolerance) + 1

                    best_match_dist = max(conserved_residue.n_term_tolerance, conserved_residue.c_term_tolerance) + 1
                    match_index = -1
                    for j in range(min_index, max_index+1):
                        dist = abs(i - j)
                        if alignment[1].seq[j] in conserved_residue.allowed_aa and dist < best_match_dist:
                            best_match_dist = dist
                            match_index = candidate_seq_index
                            last_conserved_index_match = i

                    if match_index == -1:
                        return

        # if the candidate passes all criteria, write the alignment file to disk
        with open(aln_path, 'w') as aln_file_permanent:
            aln_file.seek(0)
            shutil.copyfileobj(aln_file, aln_file_permanent)


########################################################################################################################
# Parse command line arguments
########################################################################################################################

# input files
seed_fasta_dir = sys.argv[1]
seed_index_file = open(sys.argv[2])
db_fasta_file = open(sys.argv[3])

# alignment parameters
min_alignment_identity_threshold = float(sys.argv[4])

# job parameters
job_id = int(sys.argv[5]) - 1
num_jobs = int(sys.argv[6])

# output
out_path = sys.argv[7]
index_outfile = open(sys.argv[8], 'w')

########################################################################################################################
# Initialize
########################################################################################################################

job_out_path = os.path.join(out_path, str(job_id))

# create directory for job output
if not os.path.exists(job_out_path):
    os.mkdir(job_out_path)

# read fasta files
db_sequences = SeqIO.parse(db_fasta_file, 'fasta')

# create seed proteins
seed_id2proteins = dict()
# parse conserved residues
for line in seed_index_file:
    if line.startswith("PROTEIN="):
        seed_id = line.rstrip().split("=")[-1]
        seed_fasta_path = os.path.join(seed_fasta_dir, util.fasta_id_to_fasta_filename(seed_id))
        with open(seed_fasta_path) as seed_fasta_file:
            seed_seq = next(SeqIO.parse(seed_fasta_file, 'fasta'))
        seed_id2proteins[seed_id] = SeedProtein(seed_seq, seed_fasta_path)
    else:
        [seed_aa, external_index, allowed_aa, n_term_tolerance, c_term_tolerance] = line.rstrip().split(",")

        seed_id2proteins[seed_id].add_conserved_residue(seed_aa, int(external_index), allowed_aa,
                                                        int(n_term_tolerance), int(c_term_tolerance))



########################################################################################################################
# Align
########################################################################################################################

start_time = time.perf_counter()

for i, candidate in enumerate(db_sequences):
    if (i % num_jobs) == job_id:
        align_all_seeds(seed_id2proteins, candidate, out_path, job_id, min_alignment_identity_threshold)

end_time = time.perf_counter()

print(start_time, end_time, end_time - start_time)





