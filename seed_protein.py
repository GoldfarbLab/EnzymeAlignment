import sys


class ConservedResidue:
    def __init__(self, seed_aa, external_index, allowed_aa, n_term_tolerance, c_term_tolerance):
        self.seed_aa = seed_aa
        self.index = external_index - 1
        self.allowed_aa = allowed_aa
        self.n_term_tolerance = n_term_tolerance
        self.c_term_tolerance = c_term_tolerance


class MatchedResidue:
    def __init__(self, original_residue, original_index, matched_residue, matched_index, offset):
        self.original_residue = original_residue
        self.original_external_index = original_index + 1
        self.matched_residue = matched_residue
        self.matched_external_index = matched_index + 1
        self.offset = offset


class SeedProtein:
    def __init__(self, seed, fasta_path):
        self.id = seed.id
        self.seq = seed.seq
        self.fasta_path = fasta_path
        self.index2conserved_residue = dict()

    def add_conserved_residue(self, seed_aa, external_index, allowed_aa, n_term_tolerance, c_term_tolerance):
        index = external_index - 1
        if self.seq[index] != seed_aa:
            raise ValueError('Conserved residue is incorrect for {seed_id}.\n'
                             'Excepted {expected}, observed {observed}'
                             .format(seed_id=self.id,
                                     expected=seed_aa + str(index),
                                     observed=self.seq[index] + str(external_index)))

        if len(self.index2conserved_residue) > 0 and index <= max(self.index2conserved_residue):
            raise ValueError('Conserved residue should be in increasing order.\n'
                             '{seed_id}: index {current_index} came after {max_index}'
                             .format(seed_id=self.id,
                                     current_index=index + 1,
                                     max_index=max(self.index2conserved_residue) + 1))

        self.index2conserved_residue[index] = ConservedResidue(seed_aa, index, allowed_aa, n_term_tolerance, c_term_tolerance)


