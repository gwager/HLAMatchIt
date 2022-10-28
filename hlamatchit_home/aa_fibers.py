#!/usr/bin/env python
#
#
#
#
# aa_matching.py - Module for amino acid matching functions

from collections import defaultdict
from os import sep
from pickle import FALSE, TRUE
from unicodedata import name
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import random
#from "./AminoAcidMatching/aa_matching_msf" import *
#aa_mm = AAMatch(dbversion=3420)



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

#race = ["AAFA","AFA", "AFB", "AINDI", "AISC","ALANAM" , "AMIND", "API", "CARB","CARHIS", "CARIBI", "CAU", "EURCAU", "FILII","HAWI","HIS", "JAPI", "KORI", "MENAFC", "MSWHIS", "NAM", "NCHI", "SCAHIS", "SCAMB", "SCSEAI", "VIET"]

#antigen_allele_list = {}
#anti_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/OPTN_antigens_to_alleles_CPRA.txt"
#antifile = open(anti_filename, "r")

anti_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/OPTN_antigens_to_alleles_CPRA.txt"
antifile = open(anti_filename, "r")

def freqfileselect(race):
	freqs_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/freqs_9loc_2020_trim/freqs.", race,".csv"
	#freqsfile = open(freqs_filename, "r")
	return freqs_filename

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

#identifying and pulling alleles from antigens

# get the AA mature protein subsequence of any HLA allele
# Python strings start at zero, but amino acid sequences start at position 1
def getAAsubstring(allele,start_position,end_position):
	start_position = int(start_position)
	end_position = int(end_position)
	string = HLA_seq[allele].seq[start_position-1:end_position]
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
		return int(True), AA_allele1, AA_allele2 
	else:
		return int(False), AA_allele1, AA_allele2

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


#new function to get AA string mismatches from alleles
def getAAstringmatch(allele1, allele2, locus):
	start_position = ard_start_pos[locus]
	end_position = ard_end_pos[locus]
	string1 = HLA_seq[allele1].seq[start_position-1:end_position]
	string2 = HLA_seq[allele2].seq[start_position-1:end_position]
	mm_count = 0
	pos_list = []
	pos1_list = []
	pos2_list = []
	string1 = str(string1)
	string2 = str(string2)
	start_position = int(start_position)
	end_position = int(end_position)
	for pos in range(start_position,end_position):
		aa1 = string1[pos-1]
		aa2 = string2[pos-1]
		#print(aa1, '+', aa2)
		if(aa1 == aa2):
			continue
		else:
			mm_count+=1
			pos_list.append(pos)
			pos1_list.append(pos)
			pos2_list.append(pos)
			pos1_list.append(aa1)
			pos2_list.append(aa2)

	string1 = str(string1)
	string2 = str(string2)
	pos_list = ', '.join(str(item) for item in pos_list)
	pos1_list = ', '.join(str(item) for item in pos1_list)
	pos2_list = ', '.join(str(item) for item in pos2_list)
	start_position = int(start_position)
	end_position = int(end_position)
	return string1, string2, mm_count, pos_list, pos1_list, pos2_list, end_position

#definitions to generate desired webtool

def highfreq(race, allele):
	freqs_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/freqs_9loc_2020_trim/freqs." +race + ".csv"
	#print(freqs_filename)
	freqsfile = open(freqs_filename, "r")
	possible_freqs = {}
	possible_string_value = 0
	for line in freqsfile:
		#print(line)
		(Haplo,Count,Freq) = line.split(",")
		if (allele in Haplo):
			#print(allele)
			(A,C,B,DRB345,DRB1,DQA1,DQB1,DPA1,DPB1) = Haplo.split ('~')
			freq= float(Freq)
			possible_string_value += freq

	return possible_string_value

#utilized code from https://github.com/lgragert/unos-cpra-calculator/blob/main/unos/cpra/cpra2022_nonARD_lib.py
def antigen2allele(antigen):
	#print(antigen)
	anti_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/OPTN_antigens_to_alleles_CPRA.txt"
	antifile = open(anti_filename, "r")
	antigen_allele_list = {}
	for row in antifile:
		row=row.rstrip()
		#print(antigen)
		if antigen in row:
			(antigen,alleles_gl) = row.split("\t")
			# change list of alleles in GL string to comma-delimited
			allele_list = alleles_gl.split("/")
			antigen_allele_list[antigen] = ",".join(allele_list)
		else:
			continue
	#print(antigen)
	#print(allele_list)
	alleles = allele_list

	#print(alleles)
	return alleles

def antigen2HFallele(race,antigen):
	alleles = antigen2allele(antigen)
	#print(alleles)
	possible_alleles = defaultdict(dict)
	for allele in alleles:
		#print(allele)
		possible = highfreq(race,allele)
		if allele not in possible_alleles:
			possible_alleles[allele]=possible
		else:
			continue
	#print("antigen2alles:", possible_alleles)
	max_freq = max(possible_alleles.values())
	max_allele = max(possible_alleles, key=possible_alleles.get)
	#print(max_allele, max_freq)
	return max_allele, max_freq

def fibers(locus,pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.split(',')
	hazard = 0
	for pos in pos_list:
		print(pos)
		if (locus == "A"):
			if (pos == 12|44|63|105|111|114|152|161|166|167):
				hazard = 1.09
		if (locus== "C"):
			if (pos== 11|35|30|31|39|41|46|65|70|97|108|122|143|156|160|163|176|179):
				hazard = 1.04
		if (locus == "B"):
			if (pos == 23|24|46|67|136|145):
				hazard = 1.04
		if (locus == "DQB1"):
			if (pos == 18|49|55|66|74):
				hazard = 1.07
		if (locus == "DRB1"):
			print(locus)
			if (pos== "11") or (pos == "14") or (pos=="16") or (pos=="23") or (pos=="26") or (pos=="28") or (pos=="30") or (pos=="32") or (pos=="37") or (pos=="50") or (pos=="51") or (pos=="60") or (pos== "78"):
				hazard = 1.11
				print(hazard)
	
	return hazard




# weighted choice from https://scaron.info/blog/python-weighted-choice.html
def weighted_choice(seq, weights):
    assert len(weights) == len(seq)
    assert abs(1. - sum(weights)) < 1e-6

    x = random.random()
    for i, elmt in enumerate(seq):
        if x <= weights[i]:
            return elmt
        x -= weights[i]
