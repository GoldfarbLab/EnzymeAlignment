#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO, AlignIO
from math import floor, ceil
from seed_protein import SeedProtein, MatchedResidue
import util
import time
import tempfile
import shutil


########################################################################################################################
# Functions
########################################################################################################################

def align_all_seeds(seed_id2proteins, candidate, min_alignment_identity_threshold,
                    alignment_root_out_path, residues_root_out_path, fasta_root_out_path):

    seed_id2residues = dict()
    processed_seed_ids = set()
    seed_ids_grouped_by_residues = []

    for seed in seed_id2proteins.values():
        if (len(candidate.seq) * min_alignment_identity_threshold > len(seed.seq) or  # too long
                len(candidate.seq) < len(seed.seq) * min_alignment_identity_threshold):  # too short
            continue

        matching_residues = align_to_seed(seed, candidate,
                                          alignment_root_out_path, min_alignment_identity_threshold)

        if matching_residues:
            seed_id2residues[seed.id] = matching_residues

    # group seeds together if they led to identical conserved residues. This cleans up the output.
    for seed_id in seed_id2residues:
        if seed_id in processed_seed_ids:
            continue

        seed_ids_grouped_by_residues.append([[seed_id], seed_id2residues[seed_id]])
        processed_seed_ids.add(seed_id)

        for seed_id2 in seed_id2residues:
            if seed_id2 in processed_seed_ids:
                continue

            if seed_id2residues[seed_id2] == seed_id2residues[seed_id]:
                seed_ids_grouped_by_residues[-1][0].append(seed_id2)
                processed_seed_ids.add(seed_id2)

    # write results: fasta and conserved residues
    if len(seed_id2residues) > 0:
        SeqIO.write(candidate, os.path.join(fasta_root_out_path, candidate.id + ".fasta"), 'fasta')
        write_matching_conserved_residues(candidate, seed_ids_grouped_by_residues,
                                          os.path.join(residues_root_out_path, candidate.id + "_conserved_residues.txt"))


def align_to_seed(seed, candidate, alignment_root_out_path, min_alignment_identity_threshold):

    aln_path = os.path.join(alignment_root_out_path,
                            util.fasta_id_to_filename(candidate.id) + "_vs_" + util.fasta_id_to_filename(seed.id) + ".aln.txt")

    with tempfile.SpooledTemporaryFile(max_size=1e6, mode="w+") as aln_file:
        p = subprocess.Popen(["stretcher",
                              "-bsequence", seed.fasta_path,
                              "-gapopen", "1",
                              "-gapextend", "1",
                              "-aformat3", "srspair",
                              "-auto",
                              "-filter",
                              "-sprotein1",
                              "-sprotein2"], stdout=aln_file, stdin=subprocess.PIPE, universal_newlines=True)

        p.communicate(input=util.protein2fasta_str(candidate))

        # check if the candidate passes the sequence identity threshold
        aln_file.seek(0)
        for line in aln_file:
            if line.startswith("# Identity: "):
                identity = float(line.rstrip().split("(")[-1][0:-2]) / 100.0
                if identity < min_alignment_identity_threshold:
                    return False
                break

        # check for conserved residues
        aln_file.seek(0)
        alignment = AlignIO.read(aln_file, "emboss")

        seed_seq_index = -1
        candidate_seq_index = -1
        last_conserved_index_match = -1

        matched_residues = []

        seed_aln_seq = alignment[1].seq
        candidate_aln_seq = alignment[0].seq
        aln_length = len(seed_aln_seq)

        for i in range(0, aln_length):
            if candidate_aln_seq[i] != "-":
                candidate_seq_index += 1
            if seed_aln_seq[i] != "-":
                seed_seq_index += 1
                if seed_seq_index in seed.index2conserved_residue:
                    conserved_residue = seed.index2conserved_residue[seed_seq_index]

                    # Find closest match in allowable subsequence. Tie breaks go to the left.
                    min_index = max(0, i - conserved_residue.n_term_tolerance, last_conserved_index_match)
                    max_index = min(aln_length-1, i + conserved_residue.c_term_tolerance) + 1

                    best_match_dist = max(conserved_residue.n_term_tolerance, conserved_residue.c_term_tolerance) + 1
                    match_index = -1
                    for j in range(min_index, max_index+1):
                        dist = abs(i - j)
                        if candidate_aln_seq[j] in conserved_residue.allowed_aa and dist < best_match_dist:
                            best_match_dist = dist
                            match_index = j
                            last_conserved_index_match = i

                    if match_index == -1:
                        return False
                    else:
                        matched_residues.append(MatchedResidue(conserved_residue.seed_aa, conserved_residue.index,
                                                               candidate_aln_seq[match_index], match_index,
                                                               match_index - i))

        # if the candidate passes all criteria, write the alignment file to disk
        with open(aln_path, 'w') as aln_file_permanent:
            aln_file.seek(0)
            shutil.copyfileobj(aln_file, aln_file_permanent)

            return matched_residues


def write_matching_conserved_residues(candidate, seed_ids_grouped_by_residues, out_path):
    out_file = open(out_path, "w")

    out_file.write("CANDIDATE: " + candidate.id + "\n")

    for group in seed_ids_grouped_by_residues:
        out_file.write("SEEDS: " + ", ".join(group[0]) + "\n")

        for cr in group[1]:
            out_file.write(cr.original_residue + str(cr.original_external_index) + " -> " +
                           cr.matched_residue + str(cr.matched_external_index) + " OFFSET: " + str(cr.offset) + "\n")


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
alignment_root_out_path = sys.argv[7]
residues_root_out_path = sys.argv[8]
fasta_root_out_path = sys.argv[9]


########################################################################################################################
# Initialize
########################################################################################################################

# create directory for job output
if not os.path.exists(alignment_root_out_path):
    os.mkdir(alignment_root_out_path)

if not os.path.exists(residues_root_out_path):
    os.mkdir(residues_root_out_path)

if not os.path.exists(fasta_root_out_path):
    os.mkdir(fasta_root_out_path)

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
        align_all_seeds(seed_id2proteins, candidate, min_alignment_identity_threshold,
                        alignment_root_out_path, residues_root_out_path, fasta_root_out_path)

end_time = time.perf_counter()

print(start_time, end_time, end_time - start_time)





