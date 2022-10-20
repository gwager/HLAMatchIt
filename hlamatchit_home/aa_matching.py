#!/usr/bin/env python

# aa_matching.py - Module for amino acid matching functions

from pickle import FALSE, TRUE
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import random

# mature protein sequence offsets differ by locus
# get offsets by examining mature protein length in IMGT/HLA alignment tool
# https://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/align.cgi - Mature protein
hlaProteinOffset = {
    "A" : 24, # 365 versus 341 mature
    "B" : 24,
    "C" : 24,
    "DRA" : 25,
    "DRB1" : 29,
    "DRB3" : 29,
    "DRB4" : 29,
    "DRB5" : 29,
    "DQA1" : 23,
    "DQB1" : 32,
    "DPA1" : 31,
    "DPB1" : 29,
}

def getMatureProteinOffset(locus):
        return hlaProteinOffset.get(locus, "Invalid HLA Locus")

# use only ARD positions
ard_start_pos = {
    "A" : 1,
    "B" : 1,
    "C" : 1,
    "DRB1" : 1,
	"DRB3" : 1,
	"DRB4" : 1,
	"DRB5" : 1,	
    "DQA1" : 1,
    "DQB1" : 1,
    "DPA1" : 1,
    "DPB1" : 1,
}
ard_end_pos = {
    "A" : 182,
    "B" : 182,
    "C" : 182,
    "DRB1" : 94,
	"DRB3" : 94,
	"DRB4" : 94,
	"DRB5" : 94,	
    "DQA1" : 94,
    "DQB1" : 94,
    "DPA1" : 94,
    "DPB1" : 94,
}

# lots of incomplete ARD sequences in IMGT/HLA
# this should be handled in reference alignment, not here
ard_start_pos_incomplete = {
    "A" : 2, # A*02:50
    "B" : 2, # B*07:30
    "C" : 2, # C*01:10
    "DRB1" : 7, # DRB1*08:19
	"DRB3" : 2, # TBD
	"DRB4" : 2, # TBD
	"DRB5" : 2, # TBD
    "DQA1" : 6, # DQA1*01:06
    "DQB1" : 6, # DQB1*05:100
    "DPA1" : 11, # DPA1*01:03:02
    "DPB1" : 6, # DPB1*01:01:03
}
ard_end_pos_incomplete = {
    "A" : 182, # A*02:50
    "B" : 182, # B*07:30
    "C" : 182, # C*01:10
    "DRB1" : 92, # DRB1*08:19
	"DRB3" : 2, # TBD
	"DRB4" : 2, # TBD
	"DRB5" : 2, # TBD
    "DQA1" : 87, # DQA1*01:06
    "DQB1" : 94, # DQB1*05:100
    "DPA1" : 84, # DPA1*01:03:02
    "DPB1" : 92, # DPB1*01:01:03
}


# load protein sequence file to get full protein sequences into SeqRecord file
seq_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/IMGT_HLA_Full_Protein_3330.txt"
seqfile = open(seq_filename, "r")

HLA_full_allele = {} # Full four-field allele names
HLA_seq = {} # Two-field
for line in seqfile:
	(allele,full_protein) = line.split("\t")
	if allele == "Allele":
		continue

	# skip missing sequences in IMGT/HLA file
	if (len(full_protein) <10): # #NAME?
		print ("Missing Sequence:" + allele)
		continue

	(imgt_loc,full_allele) = allele.split("*")
	(gene,loc) = imgt_loc.split('-')
	mature_protein = full_protein[getMatureProteinOffset(loc):]
	# print (allele + " " + mature_protein)
	loc_full_allele = loc + "*" + full_allele

	allele_fields = full_allele.split(':')
	two_field_allele = allele_fields[0] + ":" + allele_fields[1]
	loc_two_field_allele = loc + "*" + two_field_allele

	record = SeqRecord(Seq(mature_protein,IUPAC.protein), mature_protein, '','')
	# print (record)

	# full allele name
	HLA_full_allele[loc_full_allele] = record

	# don't overwrite two-field alleles with new sequences - more likely to be incomplete
	if (loc_two_field_allele not in HLA_seq):
		HLA_seq[loc_two_field_allele] = record

	# print (HLA_seqrecord_dict[allele])
	
	# TODO - add feature annotation to SeqIO object
	# https://biopython.org/wiki/SeqRecord
	# e.g. - which allele contain Bw4/Bw6 epitopes, bind LILRB1, etc
	# https://www.biostars.org/p/57549/

	# print (HLA_seqrecord_dict[allele].seq)


# get the AA mature protein subsequence of any HLA allele
# Python strings start at zero, but amino acid sequences start at position 1
def getAAsubstring(allele,start_position,end_position):
	return HLA_seq[allele].seq[start_position-1:end_position]

# get the AA at any position of any HLA allele
def getAAposition(allele,position):
	position = int(position)
	AA = HLA_seq[allele].seq[position-1]
	return HLA_seq[allele].seq[position-1]

# get the SFVT epitope name (with underscore)
def getEpitope(allele,position_list):
	sfvt_list = []
	for position in position_list:
		AA = getAAposition(allele,position)
		sfvt_aa = str(position) + AA
		sfvt_list.append(sfvt_aa)

	sfvt = "_".join(sfvt_list)
	return (sfvt)

# given two alleles at same locus and a position - Yes/No if mismatched
def isPositionMismatched(allele1,allele2,position):
	AA_allele1 = getAAposition(allele1,position)
	AA_allele2 = getAAposition(allele2,position)
	position = int(position)
	if (AA_allele1 != AA_allele2):
		return int(True)
	else:
		return int(False)

# count number of mismatches at position between donor and recip
# getTruth function from runMatchMC
def count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip):
	mm_count = 0
	if ((aa1_donor != aa1_recip) & (aa1_donor != aa2_recip)):
		mm_count+=1
	if ((aa2_donor != aa1_recip) & (aa2_donor != aa2_recip)):
		mm_count+=1
	return mm_count

# count number of mismatches between alleles at a given position, adjusting for donor homozygosity
def count_AA_Mismatches_Allele(allele1_donor,allele2_donor,allele1_recip,allele2_recip,position):
	donor_homoz = 0
	if (allele1_donor == allele2_donor):
		donor_homoz = 1
	aa1_donor = getAAposition(allele1_donor,position)
	aa2_donor = getAAposition(allele2_donor,position)
	aa1_recip = getAAposition(allele1_recip,position)
	aa2_recip = getAAposition(allele2_recip,position)

	mm_count = count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

	if ((mm_count == 2) & (donor_homoz == 1)):
		mm_count = 1

	return mm_count

def count_AA_Mismatches_SFVT(allele1_donor,allele2_donor,allele1_recip,allele2_recip,position_list):
	donor_homoz = 0
	mm_total = 0
	if (allele1_donor == allele2_donor):
		donor_homoz = 1
	
	# increment MM count for all positions in list
	for position in position_list:
		aa1_donor = getAAposition(allele1_donor,position)
		aa2_donor = getAAposition(allele2_donor,position)
		aa1_recip = getAAposition(allele1_recip,position)
		aa2_recip = getAAposition(allele2_recip,position)

		mm_count = count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

		if ((mm_count == 2) & (donor_homoz == 1)):
			mm_count = 1

		mm_total = mm_total + mm_count

	return mm_total


def AA_MM(aa1_donor,aa2_donor,aa1_recip,aa2_recip):
	mm_count = 0
	if ((aa1_donor != aa1_recip) & (aa1_donor != aa2_recip)):
		mm_count=1
	if ((aa2_donor != aa1_recip) & (aa2_donor != aa2_recip)):
		mm_count=1
	return mm_count


# any there any mismatches between alleles at a given position
def AA_MM_Allele(allele1_donor,allele2_donor,allele1_recip,allele2_recip,position):
	aa1_donor = getAAposition(allele1_donor,position)
	aa2_donor = getAAposition(allele2_donor,position)
	aa1_recip = getAAposition(allele1_recip,position)
	aa2_recip = getAAposition(allele2_recip,position)

	is_mm = AA_MM(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

	return is_mm

# weighted choice from https://scaron.info/blog/python-weighted-choice.html
def weighted_choice(seq, weights):
    assert len(weights) == len(seq)
    assert abs(1. - sum(weights)) < 1e-6

    x = random.random()
    for i, elmt in enumerate(seq):
        if x <= weights[i]:
            return elmt
        x -= weights[i]
