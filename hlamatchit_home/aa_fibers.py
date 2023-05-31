#!/usr/bin/env python
#
#
#
#
# aa_matching.py - Module for amino acid matching functions

from collections import defaultdict
from itertools import islice
from itertools import permutations
from itertools import product
from os import sep
from pickle import FALSE, TRUE
from unicodedata import name
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import random
from .aa_matching_msf import *
import timeit
import gzip
aa_mm = AAMatch(dbversion=3500)
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import psutil
import collections
import numbers
#Actual binned pos
#dr_fibers_pos =  [10,11,12,13,14,26,32,33,37,73,77]
#dq_fibers_pos = [30,53,55,71,74,77,89]
#mem = psutil.virtual_memory()
#top 4 with additional high ald pos
dr_fibers_pos =  [13,26]
dq_fibers_pos = [30,55]

#Create antigen dict from OPTN
anti_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/OPTN_antigens_to_alleles_CPRA.txt"
#antifile = open(anti_filename, "r")
antigen_dict = defaultdict(list)
with open(anti_filename, 'r') as f:
	for line in f:
		(Antigen, Alleles)= line.rstrip("\n").split('\t')
		if (Antigen == "OPTN_Antigen"):
			continue
		allele_string = Alleles.split("/")
		#alleles = ",".join(allele_string)
		#allele_list = list(allele_string.split(","))
		#allele_list = allele_string.split(",").split(",")
		#print(Antigen)
		#print(allele_string)
		#for allele in allele_list:
			#antigen_dict[Antigen]+allele
		antigen_dict[Antigen].extend(allele_string)
		#d1=dict.fromkeys(allele_string,Antigen)
#antigen_dict[0]
#print(antigen_dict)
	#for row in f:
		#row=row.rstrip()
		#print(antigen)
		#if antigen in row:
		#(antigen,alleles_gl) = row.split("\t")
		# change list of alleles in GL string to comma-delimited
		#allele_string = alleles_gl.split("/")
		#alleles = ",".join(allele_string)
		#allele_list = alleles.split(",")
		#print(allele_list)
 
		#(Antigen, Alleles) = line.split('\t')
		#if (Antigen == "OPTN_Antigen"):
			#continue
		#items = line.split()
		#key, values = Antigen, Alleles
		#allele_list = Alleles.split('/')
		
print("Antigen Dictionary Generated")
#print(antigen_dict)

	
#print(freqsfile)
			#antigen_allele_list[antigen] = ",".join(allele_list)
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

#Generte and Store all string from DR and DQ
drloc = "DRB1"
dqloc = "DQB1"
dr_start_position = ard_start_pos[drloc]
dr_end_position = ard_end_pos[drloc]

dq_start_position = ard_start_pos[dqloc]
dq_end_position = ard_end_pos[dqloc]

#allele_list = allele_string.split(",")
drb_string_dict = defaultdict(dict)
dqb_string_dict = defaultdict(dict)

drb_list = [ v for k, v in antigen_dict.items() if 'DR' in k]
#drb_alleles = drb_list.split("]")
#print(drb_list)
for drb_allele_list in drb_list:
	for allele in drb_allele_list:
		if (allele == "DRB3*03:22"):
			continue
		else:
	#print(allele)
			drb_allele_string= aa_mm.HLA_seq[allele].seq[dr_start_position-1:dr_end_position]
			drb_string_dict[allele]=drb_allele_string
print("DRB String Dict Generated")

dqb_list = [ v for k, v in antigen_dict.items() if 'DQ' in k]

for dqb_allele_list in dqb_list:
	for allele in dqb_allele_list:
		if (allele == "DRB3*03:22"):
			continue
		else:
	#print(allele)
			dqb_allele_string= aa_mm.HLA_seq[allele].seq[dq_start_position-1:dq_end_position]
			dqb_string_dict[allele]=dqb_allele_string
print("DQB String Dict Generated")
#print("DQB String Dict Generated")
#print(dict(islice(dqb_string_dict.items(), 0, 4)))

#anti_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/OPTN_antigens_to_alleles_CPRA.txt"
#antifile = open(anti_filename, "r")

def freqfileselect(race):
	freqs_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/freqs_9loc_2020_trim/freqs.", race,".csv"
	#freqsfile = open(freqs_filename, "r")
	return freqs_filename

# Define string sort to keep allele names as is when sorting
def sort_string(mylist):
    #list to store our result
    ls = []
    while mylist:
        #append minimum element of the list of strings
        ls.append(min(mylist))
        #remove that minimum element from original string
        mylist.remove(min(mylist))
    return ls

#Define Cumulative Distribution Function for dictionaries    
def cumulative(function_dict):
	function_dict = dict( sorted(function_dict.items(), key=lambda kv: kv[1], reverse=False))
	sum_values = sum(function_dict.values())
	#print(sum_values)
	F= defaultdict(dict)
	temp = 0
	for i in function_dict:
		j = function_dict.get(i)
		temp =  temp+j
		cdf = temp/sum_values
		#cdf = k - avg
		if (cdf >= .00001):
			if i not in F:
				F[i]=cdf
			else:
				continue
		else:
			continue
	return F

def mm_round(test_dict):
	res = defaultdict(dict)
	for key in test_dict:
		res[key] = round(test_dict[key], 3)
	return res

#Generate function to round webtool output
#https://www.learnbyexample.org/python-nested-dictionary/
def nested_mm_round(mm_dict):
	mm_round_out = defaultdict(dict)
	for pos, aas in mm_dict.items():
		#print(pos)
		#print(aas)
		for value in aas:
			prob = aas[value]
			#print(prob)
			prob = round(prob,2)
			if value not in mm_round_out[pos]:
				mm_round_out[pos][value]= prob
			else:
				continue
	return mm_round_out

# load protein sequence file to get full protein sequences into SeqRecord file
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
	string = aa_mm.HLA_seq[allele].seq[start_position-1:end_position]
	return aa_mm.HLA_seq[allele].seq[start_position-1:end_position]

# get the AA at any position of any HLA allele
def getAAposition(allele,position):
	position = int(position)
	AA = aa_mm.HLA_seq[allele].seq[position-1]
	return aa_mm.HLA_seq[allele].seq[position-1]

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
	string1 = aa_mm.HLA_seq[allele1].seq[start_position-1:end_position]
	string2 = aa_mm.HLA_seq[allele2].seq[start_position-1:end_position]
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

def getAAgenostringmatch(da1, da2, ra1, ra2, locus):
	start_position = ard_start_pos[locus]
	end_position = ard_end_pos[locus]
	dstring1 = aa_mm.HLA_seq[da1].seq[start_position-1:end_position]
	dstring2 = aa_mm.HLA_seq[da2].seq[start_position-1:end_position]
	rstring1 = aa_mm.HLA_seq[ra1].seq[start_position-1:end_position]
	rstring2 = aa_mm.HLA_seq[ra2].seq[start_position-1:end_position]
	mm_count = 0
	pos_list = []
	pos1_list = []
	dstring1 = str(dstring1)
	dstring2 = str(dstring2)
	rstring1 = str(rstring1)
	rstring2 = str(rstring2)
	start_position = int(start_position)
	end_position = int(end_position)
	for pos in range(start_position,end_position):
		daa1 = dstring1[pos-1]
		daa2 = dstring2[pos-1]
		raa1 = rstring1[pos-1]
		raa2 = rstring2[pos-1]
		daa = daa1 + daa2
		#print(daa)
		daa = (daa)
		daa = ''.join(sorted(daa, key=str.lower))
		#print(daa)
		raa = raa1 + raa2
		raa = ''.join(sorted(raa, key=str.lower))
		#print(raa , daa)
		if(raa == daa):
			continue
		else:
			mm_count+=1
			pos_list.append(pos)
			pos1_list.append(pos)
			combo = daa + "|" + raa
			pos1_list.append(combo)

	dstring1 = str(dstring1)
	dstring2 = str(dstring2)
	rstring1 = str(rstring1)
	rstring2 = str(rstring2)
	pos_list = ', '.join(str(item) for item in pos_list)
	pos1_list = ', '.join(str(item) for item in pos1_list)
	start_position = int(start_position)
	end_position = int(end_position)
	return dstring1, dstring2, rstring1, rstring2, mm_count, pos_list, pos1_list, end_position



def highfreq(race, allele):
	freqs_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/freqs_9loc_2020_trim/freqs." +race + ".csv"
	#print(freqs_filename)
	freqsfile = open(freqs_filename, "r")
	#possible_freqs = {}
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
	possible_alleles = defaultdict(dict)
	for allele in alleles:
		#print(allele)
		possible = highfreq(race,allele)
		if allele not in possible_alleles:
			possible_alleles[allele]=possible
		else:
			continue
	#print("Possible Alleles with Population Frequencies from NMDP selected by highfreq function:", possible_alleles)
	max_freq = max(possible_alleles.values())
	max_allele = max(possible_alleles, key=possible_alleles.get)
	alleles = {x:y for x,y in possible_alleles.items() if y!=0}
	nalleles = len(alleles)
	sumant =sum(possible_alleles.values())
	#print(allelefreqs)
	return max_allele, max_freq, nalleles, sumant



def antigen2HFalleleFIBERS(race,antigen):
	alleles = antigen2allele(antigen)
	possible_alleles = defaultdict(dict)
	for allele in alleles:
		#print(allele)
		possible = highfreq(race,allele)
		if allele not in possible_alleles:
			possible_alleles[allele]=possible
		else:
			continue
	#print("Possible Alleles with Population Frequencies from NMDP selected by highfreq function:", possible_alleles)
	max_freq = max(possible_alleles.values())
	max_allele = max(possible_alleles, key=possible_alleles.get)
	alleles = {x:y for x,y in possible_alleles.items() if y!=0}
	nalleles = len(alleles)
	sumant =sum(possible_alleles.values())
	#print(allelefreqs)
	return max_allele, max_freq, nalleles, sumant, alleles



def getAAgenostringmatchFIBERS(dalleles1, dalleles2, ralleles1, ralleles2, locus):
	probs = defaultdict(dict)
	pos_freqs = defaultdict(dict)
	pos_probs = []
	aa_dr_probs = defaultdict(dict)

	#dictionary for D or R genotype frequency for prob calculations
	geno_list_donor = defaultdict(dict)
	geno_list_recip = defaultdict(dict)

	#dictionary for D|R pair genotype frequency for prob calculations
	dr_geno = defaultdict(dict)
	start_position = ard_start_pos[locus]
	end_position = ard_end_pos[locus]
	
	#calculate donor and recip genotype frequencies and create a genotype dictionary for each
	for da1 in dalleles1:
		for da2 in dalleles2:
			da = (da1)+ '+' + (da2)
			#da_un =da1 + da2
			#da = '+'.join(sorted(da_un, key=str.lower))
			if (da1 == da2):
				df1 = dalleles1.get(da1)
				gfreq = df1*df1
				geno_list_donor[da]=gfreq
			else:
				df1 = dalleles1.get(da1)
				df2 = dalleles2.get(da2)
				gfreq = 2*df1*df2
				geno_list_donor[da]=gfreq

	for ra1 in ralleles1:
		for ra2 in ralleles2:
			ra = (ra1) + '+' + (ra2)
			#ra_un =ra1 + ra2
			#ra = '+'.join(sorted(ra_un, key=str.lower))
			if (ra1 == ra2):
				rf1 = ralleles1.get(ra1)
				gfreq = rf1*rf1
				geno_list_recip[ra]=gfreq
			else:
				rf1 = ralleles1.get(ra1)
				rf2 = ralleles2.get(ra2)
				gfreq = 2*rf1*rf2
				geno_list_recip[ra]=gfreq

	#Get strings for donor and recip geno types, calculate pair frequency, and store for probability calculation
	for dgeno in geno_list_donor:
		da1,da2 = dgeno.split('+')
		dfgeno = geno_list_donor.get(dgeno)
		dstring1 = aa_mm.HLA_seq[da1].seq[start_position-1:end_position]
		dstring2 = aa_mm.HLA_seq[da2].seq[start_position-1:end_position]
		dstring1 = str(dstring1)
		dstring2 = str(dstring2)
		for rgeno in geno_list_recip:
			ra1,ra2 = rgeno.split('+')
			rfgeno = geno_list_recip.get(rgeno)
			rstring1 = aa_mm.HLA_seq[ra1].seq[start_position-1:end_position]
			rstring2 = aa_mm.HLA_seq[ra2].seq[start_position-1:end_position]
			combo = (dgeno) + '|' + (rgeno)
			rstring1 = str(rstring1)
			rstring2 = str(rstring2)

			if (dgeno == rgeno):
				drfreq = rfgeno*dfgeno
				dr_geno[combo] = drfreq
			else:
				drfreq = 2*rfgeno*dfgeno
				dr_geno[combo] = drfreq

			#based on locus check for MM at select positions and if there are store combo and freq in MM dictionary for probability calcs
			if (locus == "DRB1"):
				for pos in list(dr_fibers_pos):
					daa1 = dstring1[pos-1]
					daa2 = dstring2[pos-1]
					raa1 = rstring1[pos-1]
					raa2 = rstring2[pos-1]
					daa = daa1 + daa2
					#print(daa)
					daa = ''.join(sorted(daa, key=str.lower))
					#print(daa)
					raa = raa1 + raa2
					raa = ''.join(sorted(raa, key=str.lower))
					#print(raa , daa)
						#print(geno_list)
					show = raa + '|'+ daa
					if(raa == daa):
						#print("yes:", raa,daa)
						continue
					else:
						draafreq = (rfgeno*rfgeno)+(dfgeno*dfgeno)
						if pos not in pos_freqs:
							pos_freqs[pos]= draafreq
						else: 
							pos_freqs[pos]= pos_freqs[pos]+draafreq
						if combo not in probs:
							probs[combo]=drfreq
						else:
							continue

			if (locus == "DQB1"):
				for pos in list(dq_fibers_pos):
					daa1 = dstring1[pos-1]
					daa2 = dstring2[pos-1]
					raa1 = rstring1[pos-1]
					raa2 = rstring2[pos-1]
					daa = daa1 + daa2
					#print(daa)
					daa = ''.join(sorted(daa, key=str.lower))
					#print(daa)
					raa = raa1 + raa2
					raa = ''.join(sorted(raa, key=str.lower))
					#print(raa , daa)
						#print(geno_list)
					if(raa == daa):
						#print("yes:", raa,daa)
						continue
					else:
						draafreq = (rfgeno*rfgeno)+(dfgeno*dfgeno)
						if pos not in pos_freqs:
							pos_freqs[pos]= draafreq
						else: 
							pos_freqs[pos]= pos_freqs[pos]+draafreq
						if combo not in probs:
							probs[combo]=drfreq
						else:
							continue
				

	#get probability of MM from summation of MM dict diveded by summation of all genotype freqs dict
	sumgeno = sum(dr_geno.values())
	sumpossible =sum(probs.values())
	mm_prob = sumpossible/sumgeno

	for pos in pos_freqs:
		sumpos = pos_freqs.get(pos)
		summm =sum(pos_freqs.values())
		pos_prob = sumpos/summm
		spos =str(pos) 
		spos_prob = str(pos_prob)
		pos_combo = spos+ ':' + spos_prob
		pos_probs.append(pos_combo)
	pos_probs = ', '.join(str(item) for item in pos_probs)

	return mm_prob, pos_probs



#functions for enumarating probs as shown in DRDQFIBERS.ipynb
def antigen2HFalleleEwoP(race,antigen):
	#print("Memory when initiating function:", psutil.virtual_memory())
	#start=timeit.default_timer()
	#start_all=timeit.default_timer()
	freqs_filename = "/Users/gracelord/dev/hlamatchit/hlamatchit_home/allele_freqs/classII_freqs_" +race + ".csv"
	freqsfile = open(freqs_filename, "r")
	allele_list = antigen_dict.get(antigen)
	#allele_list = allele_string.split(",")
	#print(allele_list)
	#print("Time to pull antigen from dict:",timeit.default_timer()-start)
	possible_alleles = defaultdict(dict)
	for line in freqsfile:
		(Pop,Allele,Freq) = line.split(",")
		if(Pop=="Pop"):
			continue
		for allele in allele_list:
			freq =float(Freq)
			if (allele == Allele):
				possible_alleles[allele]=freq
			else:
				continue

	#Delete alleles with zero frequencies
	alleles = {x:y for x,y in possible_alleles.items() if y!=0}
	#print("Memory from Assigning Allele Frequencies:", psutil.virtual_memory())
	#print("Seconds to select alleles:",timeit.default_timer()-start_all)
	return alleles

def DRDQFIBERSE(ddr1alleles,ddr2alleles, rdr1alleles, rdr2alleles, ddq1alleles,ddq2alleles,rdq1alleles,rdq2alleles):
	#Dictionaries to store D or R Haplotype Frequencies

	#Generate dicts to calculate allele freqs to probs and use probs throughout
	ddr1_probs = defaultdict(dict)
	ddr2_probs = defaultdict(dict)
	ddq1_probs = defaultdict(dict)
	ddq2_probs = defaultdict(dict)

	rdr1_probs = defaultdict(dict)
	rdr2_probs = defaultdict(dict)
	rdq1_probs = defaultdict(dict)
	rdq2_probs = defaultdict(dict)


	ddr1_freqs = defaultdict(dict)
	ddr2_freqs = defaultdict(dict)
	ddq1_freqs = defaultdict(dict)
	ddq2_freqs = defaultdict(dict)

	rdr1_freqs = defaultdict(dict)
	rdr2_freqs = defaultdict(dict)
	rdq1_freqs = defaultdict(dict)
	rdq2_freqs = defaultdict(dict)


	#individual dicts for single locus and multi loci D and R cals
	dr_list_donor = defaultdict(dict)
	dq_list_donor = defaultdict(dict)
	geno_list_donor = defaultdict(dict)

	dr_list_recip = defaultdict(dict)
	dq_list_recip = defaultdict(dict)
	geno_list_recip = defaultdict(dict)

	#Dictionary to store D|R Haplotype pair frequency
	db_geno = defaultdict(dict)


	#Dictionaries to store AA-MM frequencies not MM or MM
	dr_mm_freqs = defaultdict(dict)
	dq_mm_freqs = defaultdict(dict)

	#this will give haplotype D|R pair AA-MM prob
	mm_freqs = defaultdict(dict)



	#calculate donor and recip DR and DQ frequencies and then create a genotype dictionary for each
	#calculate DR loci freq from allele level freqs
	#sort dictionary ascending and apply CDF function to cut low probabilities for optomization

	#turn all allele frequencies into probabilities to use probabilities through out
	total_recip_dr1_freq = sum(rdr1alleles.values())
	total_recip_dr2_freq = sum(rdr2alleles.values())

	total_donor_dr1_freq = sum(ddr1alleles.values())
	total_donor_dr2_freq = sum(ddr2alleles.values())

	total_recip_dq1_freq = sum(rdq1alleles.values())
	total_recip_dq2_freq = sum(rdq2alleles.values())

	total_donor_dq1_freq = sum(ddq1alleles.values())
	total_donor_dq2_freq = sum(ddq2alleles.values())

	for key, value in ddr1alleles.items():
		dastring1 = drb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DRB1"
		for pos in dr_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		#put in prob dict for cdf
		if daa1 not in ddr1_probs:
			ddr1_probs[daa1] = (value/total_donor_dr1_freq)
		else:
			ddr1_probs[daa1] = ddr1_probs[daa1]+(value/total_donor_dr1_freq)
		#put in freq dict for HW
		if daa1 not in ddr1_freqs:
			ddr1_freqs[daa1] = value
		else:
			ddr1_freqs[daa1] = ddr1_freqs[daa1]+value


	for key, value in ddr2alleles.items():
		dastring1 = drb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DRB1"
		for pos in dr_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in ddr2_probs:
			ddr2_probs[daa1] = (value/total_donor_dr2_freq)
		else:
			ddr2_probs[daa1] = ddr2_probs[daa1]+(value/total_donor_dr2_freq)
		if daa1 not in ddr2_freqs:
			ddr2_freqs[daa1] = value
		else:
			ddr2_freqs[daa1] = ddr2_freqs[daa1]+value

	for key, value in ddq1alleles.items():
		dastring1 = dqb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DQB1"
		for pos in dq_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in ddq1_probs:
			ddq1_probs[daa1] = (value/total_donor_dq1_freq)
		else:
			ddq1_probs[daa1] = ddq1_probs[daa1]+(value/total_donor_dq1_freq)
		if daa1 not in ddq1_freqs:
			ddq1_freqs[daa1] = value
		else:
			ddq1_freqs[daa1] = ddq1_freqs[daa1]+value
	
	for key, value in ddq2alleles.items():
		dastring1 = dqb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DQB1"
		for pos in dq_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in ddq2_probs:
			ddq2_probs[daa1] = (value/total_donor_dq2_freq)
		else:
			ddq2_probs[daa1] = ddq2_probs[daa1]+(value/total_donor_dq2_freq)
		if daa1 not in ddq2_freqs:
			ddq2_freqs[daa1] = value
		else:
			ddq2_freqs[daa1] = ddq2_freqs[daa1]+value

	for key, value in rdr1alleles.items():
		dastring1 = drb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DRB1"
		for pos in dr_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in rdr1_probs:
			rdr1_probs[daa1] = (value/total_recip_dr1_freq)
		else:
			rdr1_probs[daa1] = rdr1_probs[daa1]+(value/total_recip_dr1_freq)
		if daa1 not in rdr1_freqs:
			rdr1_freqs[daa1] = value
		else:
			rdr1_freqs[daa1] = rdr1_freqs[daa1]+value

	for key, value in rdr2alleles.items():
		dastring1 = drb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DRB1"
		for pos in dr_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in rdr2_probs:
			rdr2_probs[daa1] = (value/total_recip_dr2_freq)
		else:
			rdr2_probs[daa1] = rdr2_probs[daa1]+(value/total_recip_dr2_freq)
		if daa1 not in rdr2_freqs:
			rdr2_freqs[daa1] = value
		else:
			rdr2_freqs[daa1] = rdr2_freqs[daa1]+value

	for key, value in rdq1alleles.items():
		dastring1 = dqb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DQB1"
		for pos in dq_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in rdq1_probs:
			rdq1_probs[daa1] = (value/total_recip_dq1_freq)
		else:
			rdq1_probs[daa1] = rdq1_probs[daa1] + (value/total_recip_dq1_freq)
		if daa1 not in rdq1_freqs:
			rdq1_freqs[daa1] = value
		else:
			rdq1_freqs[daa1] = rdq1_freqs[daa1]+value

	for key, value in rdq2alleles.items():
		dastring1 = dqb_string_dict.get(key)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		daa1 = "DQB1"
		for pos in dq_fibers_pos:
			aa = dastring1[pos-1]
			daa1 = daa1 +"_" + str(pos) + str(aa) 
		if daa1 not in rdq2_probs:
			rdq2_probs[daa1] = (value/total_recip_dq2_freq)
		else:
			rdq2_probs[daa1] = rdq2_probs[daa1]+(value/total_recip_dq2_freq)
		if daa1 not in rdq2_freqs:
			rdq2_freqs[daa1] = value
		else:
			rdq2_freqs[daa1] = rdq2_freqs[daa1]+value

	ddr1_cdf = cumulative(ddr1_probs)
	ddr2_cdf = cumulative(ddr2_probs)
	for ddra1 in ddr1_cdf:
		for ddra2 in ddr2_cdf:
			ddra = '+'.join(sort_string([ddra1,ddra2]))
			if (ddra1 == ddra2):
				ddrf1 = ddr1_freqs.get(ddra1)
				gdrfreq = ddrf1*ddrf1
				if ddra not in dr_list_donor:
					dr_list_donor[ddra]=gdrfreq
				else:
					dr_list_donor[ddra]=dr_list_donor[ddra] + gdrfreq
			else:
				df1 = ddr1_freqs.get(ddra1)
				df2 = ddr2_freqs.get(ddra2)
				gdrfreq = 2*df1*df2
				if ddra not in dr_list_donor:
					dr_list_donor[ddra]=gdrfreq
				else:
					dr_list_donor[ddra]=dr_list_donor[ddra] + gdrfreq

	ddq1_cdf = cumulative(ddq1_probs)
	ddq2_cdf = cumulative(ddq2_probs)

	#calculate DQ loci freq from allele level freqs
	for ddqa1 in ddq1_cdf:
		for ddqa2 in ddq2_cdf:
			ddqa = '+'.join(sort_string([ddqa1,ddqa2]))
			if (ddqa1 == ddqa2):
				ddqf1 = ddq1_freqs.get(ddqa1)
				gdqfreq = ddqf1*ddqf1
				if ddqa not in dq_list_donor:
					dq_list_donor[ddqa]=gdqfreq
				else:
					dq_list_donor[ddqa]=dq_list_donor[ddqa] + gdqfreq
			else:
				df1 = ddq1_freqs.get(ddqa1)
				df2 = ddq2_freqs.get(ddqa2)
				gdqfreq = 2*df1*df2
				if ddqa not in dq_list_donor:
					dq_list_donor[ddqa]=gdqfreq
				else:
					dq_list_donor[ddqa]=dq_list_donor[ddqa] + gdqfreq

	ddr_cdf = cumulative(dr_list_donor)
	ddq_cdf = cumulative(dq_list_donor)
	dr_list_donor = dict( sorted(dr_list_donor.items(), key=lambda kv: kv[1], reverse=True))	
	dq_list_donor = dict( sorted(dq_list_donor.items(), key=lambda kv: kv[1], reverse=True))	


	#gives all key pernumerations and multiplies values
	for key, value in dr_list_donor.items():
		for key2, value2 in dq_list_donor.items():
			geno_list_donor["{}^{}".format(key, key2)] = value * value2
	#print(geno_list_donor)
	#print("Memory after Donor calculations and pernumeration of multi-loci:", psutil.virtual_memory())
	geno_list_donor = dict( sorted(geno_list_donor.items(), key=lambda kv: kv[1], reverse=True))	


	rdr1_cdf = cumulative(rdr1_probs)
	rdr2_cdf = cumulative(rdr2_probs)
	
	#calculate DR loci freq from allele level freqs
	for rdra1 in rdr1_cdf:
		for rdra2 in rdr2_cdf:
			rdra = '+'.join(sort_string([rdra1,rdra2]))
			if (rdra1 == rdra2):
				rdrf1 = rdr1_freqs.get(rdra1)
				gdrfreq = rdrf1*rdrf1
				if rdra not in dr_list_recip:
					dr_list_recip[rdra]=gdrfreq
				else:
					dr_list_recip[rdra]=dr_list_recip[rdra] + gdrfreq
			else:
				df1 = rdr1_freqs.get(rdra1)
				df2 = rdr2_freqs.get(rdra2)
				gdrfreq = df1*df2*2
				if rdra not in dr_list_recip:
					dr_list_recip[rdra]=gdrfreq
				else:
					dr_list_recip[rdra]=dr_list_recip[rdra] + gdrfreq


	rdq1_cdf = cumulative(rdq1_probs)
	rdq2_cdf = cumulative(rdq2_probs)
	#calculate DQ loci freq from allele level freqs
	for rdqa1 in rdq1_cdf:
		for rdqa2 in rdq2_cdf:
			rdqa = '+'.join(sort_string([rdqa1,rdqa2]))
			if (rdqa1 == rdqa2):
				rdqf1 = rdq1_freqs.get(rdqa1)
				gdqfreq = rdqf1*rdqf1
				if rdqa not in dq_list_recip:
					dq_list_recip[rdqa]=gdqfreq
				else:
					dq_list_recip[rdqa]=dq_list_recip[rdqa] + gdqfreq
			else:
				df1 = rdq1_freqs.get(rdqa1)
				df2 = rdq2_freqs.get(rdqa2)
				gdqfreq = 2*df1*df2
				if rdqa not in dq_list_recip:
					dq_list_recip[rdqa]=gdqfreq
				else:
					dq_list_recip[rdqa]=dq_list_recip[rdqa] + gdqfreq
					
	#calculate DR|DQ geno freq from individual DR and DQ loci freqs level freqs

	#gives all key pernumerations and multiplies values
	rdr_cdf = cumulative(dr_list_recip)
	rdq_cdf = cumulative(dq_list_recip)
	dr_list_recip = dict( sorted(dr_list_recip.items(), key=lambda kv: kv[1], reverse=True))	
	dq_list_recip = dict( sorted(dq_list_recip.items(), key=lambda kv: kv[1], reverse=True))	


	for key, value in dr_list_recip.items():
		for key2, value2 in dq_list_recip.items():
			geno_list_recip["{}^{}".format(key, key2)] = value * value2


	####Adding to base off genotype probs instead of using frequencies as requested by Dr.Gragert
	total_donor_geno_freq = sum(geno_list_donor.values())
	total_recip_geno_freq = sum(geno_list_recip.values())

	#Get strings for donor and recip geno types, calculate pair frequency, and store for probability calculation
	#gives all key pernumerations and multiplies probs
	for key, value in geno_list_donor.items():
		for key2, value2 in geno_list_recip.items():
			db_geno["{}|{}".format(key, key2)] = (value/total_donor_geno_freq) * (value2/total_recip_geno_freq)

	db_prob= sum(db_geno.values())


	#need one loop for each pos of interest
	for combo in db_geno:
		drprob = db_geno.get(combo)
		dgeno,rgeno = combo.split('|')
		#print(dgeno)
		ddr,ddq = dgeno.split('^')
		#print(ddr)
		da1,da2 = ddr.split('+')
		#dq1,dq2 = ddq.split('+')

		rdr,rdq = rgeno.split('^')

		ra1,ra2 = rdr.split('+')
		#rq1,rq2 = rdq.split('+')
		

		# hard coding by number of pos must be changed if number of AA-MM pos are changed
		loci, r1pos1, r1pos2 = ra1.split('_')
		loci, r2pos1, r2pos2 = ra2.split('_')

		loci, d1pos1, d1pos2 = da1.split('_')
		loci, d2pos1, d2pos2 = da2.split('_')

		dpos1 = '_'.join(sort_string([d1pos1,d2pos1]))
		#dpos2 = d1pos2 + '_' + d2pos2


		rpos1 = '_'.join(sort_string([r1pos1,r2pos1])) 
		#rpos2 = r1pos2 + '_' + r2pos2

		#Needs to be changed to reflct number of FIBERS high risk pos
		#name,rdrpos1,rdrpos2 = rdq.split('+')
		#print(da1, da2, ra1, ra2)
		#print(dq1, dq2, rq1, rq2)
		if (d1pos1 == d2pos1):
			if (d1pos1 == r1pos1 or d1pos1 == r2pos1):
				continue
				#store D|R AA assignments and probability of Matches
			else:
				show = dpos1 + "|" + rpos1
				#store D|R AA assignments and probability of MisMatches
				if show not in dr_mm_freqs:
					dr_mm_freqs[show] = drprob
				else: 
					dr_mm_freqs[show] = dr_mm_freqs[show] + drprob
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#Else statement if Donor AA assignment is heterozygous
		else:
			if (d1pos1 == r1pos1 or d1pos1 == r2pos1):
				if (d2pos1 == r1pos1 or d2pos1 == r2pos1):
					continue
				else:
					show = dpos1 + "|" + rpos1
					#store D|R AA assignments and probability of Mismatches
					if show not in dr_mm_freqs:
						dr_mm_freqs[show] = drprob
					else: 
						dr_mm_freqs[show] = dr_mm_freqs[show] + drprob		
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue
			else:
				show = dpos1 + "|" + rpos1
				#store D|R AA assignments and probability of Mismatches
				if show not in dr_mm_freqs:
					dr_mm_freqs[show] = drprob
				else: 
					dr_mm_freqs[show] = dr_mm_freqs[show] + drprob		
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#print(type(dq1))

	for combo in db_geno:
		drprob = db_geno.get(combo)
		dgeno,rgeno = combo.split('|')
		#print(dgeno)
		ddr,ddq = dgeno.split('^')
		#print(ddr)
		da1,da2 = ddr.split('+')
		#dq1,dq2 = ddq.split('+')

		rdr,rdq = rgeno.split('^')

		ra1,ra2 = rdr.split('+')
		#rq1,rq2 = rdq.split('+')
		

		# hard coding by number of pos must be changed if number of AA-MM pos are changed
		loci, r1pos1, r1pos2 = ra1.split('_')
		loci, r2pos1, r2pos2 = ra2.split('_')

		loci, d1pos1, d1pos2 = da1.split('_')
		loci, d2pos1, d2pos2 = da2.split('_')

		#dpos1 = d1pos1 + '_' + d2pos1
		dpos2 = '_'.join(sort_string([d1pos2,d2pos2]))


		#rpos1 = r1pos1 + '_' + r2pos1
		rpos2 = '_'.join(sort_string([r1pos2,r2pos2]))

		#Needs to be changed to reflct number of FIBERS high risk pos
		#name,rdrpos1,rdrpos2 = rdq.split('+')
		#print(da1, da2, ra1, ra2)
		#print(dq1, dq2, rq1, rq2)
		if (d1pos2 == d2pos2):
			if (d1pos2 == r1pos2 or d1pos2 == r2pos2):
				continue
				#store D|R AA assignments and probability of Matches
			else:
				show = dpos2 + "|" + rpos2
				#store D|R AA assignments and probability of MisMatches
				if show not in dr_mm_freqs:
					dr_mm_freqs[show] = drprob
				else: 
					dr_mm_freqs[show] = dr_mm_freqs[show] + drprob
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#Else statement if Donor AA assignment is heterozygous
		else:
			if (d1pos2 == r1pos2 or d1pos1 == r2pos2):
				if (d2pos2 == r1pos2 or d2pos2 == r2pos2):
					continue
				else:
					show = dpos2 + "|" + rpos2
					#store D|R AA assignments and probability of Mismatches
					if show not in dr_mm_freqs:
						dr_mm_freqs[show] = drprob
					else: 
						dr_mm_freqs[show] = dr_mm_freqs[show] + drprob		
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue
			else:
				show = dpos2 + "|" + rpos2
				#store D|R AA assignments and probability of Mismatches
				if show not in dr_mm_freqs:
					dr_mm_freqs[show] = drprob
				else: 
					dr_mm_freqs[show] = dr_mm_freqs[show] + drprob		
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue

	#for DQ Pos
	for combo in db_geno:
		drprob = db_geno.get(combo)
		dgeno,rgeno = combo.split('|')
		#print(dgeno)
		ddr,ddq = dgeno.split('^')
		#print(ddr)
		#da1,da2 = ddr.split('+')
		dq1,dq2 = ddq.split('+')

		rdr,rdq = rgeno.split('^')

		#ra1,ra2 = rdr.split('+')
		rq1,rq2 = rdq.split('+')
		

		# hard coding by number of pos must be changed if number of AA-MM pos are changed
		loci, r1pos1, r1pos2 = rq1.split('_')
		loci, r2pos1, r2pos2 = rq2.split('_')

		loci, d1pos1, d1pos2 = dq1.split('_')
		loci, d2pos1, d2pos2 = dq2.split('_')

		dpos1 = '_'.join(sort_string([d1pos1,d2pos1]))
		#dpos2 = d1pos2 + '_' + d2pos2


		rpos1 = '_'.join(sort_string([r1pos1,r2pos1]))
		#rpos2 = r1pos2 + '_' + r2pos2

		#Needs to be changed to reflct number of FIBERS high risk pos
		#name,rdrpos1,rdrpos2 = rdq.split('+')
		#print(da1, da2, ra1, ra2)
		#print(dq1, dq2, rq1, rq2)
		if (d1pos1 == d2pos1):
			if (d1pos1 == r1pos1 or d1pos1 == r2pos1):
				continue
				#store D|R AA assignments and probability of Matches
			else:
				show = dpos1 + "|" + rpos1
				#store D|R AA assignments and probability of MisMatches
				if show not in dq_mm_freqs:
					dq_mm_freqs[show] = drprob
				else: 
					dq_mm_freqs[show] = dq_mm_freqs[show] + drprob
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#Else statement if Donor AA assignment is heterozygous
		else:
			if (d1pos1 == r1pos1 or d1pos1 == r2pos1):
				if (d2pos1 == r1pos1 or d2pos1 == r2pos1):
					continue
				else:
					show = dpos1 + "|" + rpos1
					#store D|R AA assignments and probability of Mismatches
					if show not in dq_mm_freqs:
						dq_mm_freqs[show] = drprob
					else: 
						dq_mm_freqs[show] = dq_mm_freqs[show] + drprob		
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue
			else:
				show = dpos1 + "|" + rpos1
				#store D|R AA assignments and probability of Mismatches
				if show not in dq_mm_freqs:
					dq_mm_freqs[show] = drprob
				else: 
					dq_mm_freqs[show] = dq_mm_freqs[show] + drprob		
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#print(type(dq1))

	for combo in db_geno:
		drprob = db_geno.get(combo)
		dgeno,rgeno = combo.split('|')
		#print(dgeno)
		ddr,ddq = dgeno.split('^')
		#print(ddr)
		#da1,da2 = ddr.split('+')
		dq1,dq2 = ddq.split('+')

		rdr,rdq = rgeno.split('^')

		#ra1,ra2 = rdr.split('+')
		rq1,rq2 = rdq.split('+')
		

		# hard coding by number of pos must be changed if number of AA-MM pos are changed
		loci, r1pos1, r1pos2 = rq1.split('_')
		loci, r2pos1, r2pos2 = rq2.split('_')

		loci, d1pos1, d1pos2 = dq1.split('_')
		loci, d2pos1, d2pos2 = dq2.split('_')

		#dpos1 = d1pos1 + '_' + d2pos1
		dpos2 = '_'.join(sort_string([d1pos2,d2pos2]))


		#rpos1 = r1pos1 + '_' + r2pos1
		rpos2 = '_'.join(sort_string([r1pos2,r2pos2]))

		#Needs to be changed to reflct number of FIBERS high risk pos
		#name,rdrpos1,rdrpos2 = rdq.split('+')
		#print(da1, da2, ra1, ra2)
		#print(dq1, dq2, rq1, rq2)
		if (d1pos2 == d2pos2):
			if (d1pos2 == r1pos2 or d1pos2 == r2pos2):
				continue
				#store D|R AA assignments and probability of Matches
			else:
				show = dpos2 + "|" + rpos2
				#store D|R AA assignments and probability of MisMatches
				if show not in dq_mm_freqs:
					dq_mm_freqs[show] = drprob
				else: 
					dq_mm_freqs[show] = dq_mm_freqs[show] + drprob
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
		#Else statement if Donor AA assignment is heterozygous
		else:
			if (d1pos2 == r1pos2 or d1pos2 == r2pos2):
				if (d2pos2 == r1pos2 or d2pos2 == r2pos2):
					continue
				else:
					show = dpos2 + "|" + rpos2
					#store D|R AA assignments and probability of Mismatches
					if show not in dq_mm_freqs:
						dq_mm_freqs[show] = drprob
					else: 
						dq_mm_freqs[show] = dq_mm_freqs[show] + drprob		
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue
			else:
				show = dpos2 + "|" + rpos2
				#store D|R AA assignments and probability of Mismatches
				if show not in dq_mm_freqs:
					dq_mm_freqs[show] = drprob
				else: 
					dq_mm_freqs[show] = dq_mm_freqs[show] + drprob		
				#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
				if combo not in mm_freqs:
					mm_freqs[combo] = drprob
				else:
					continue
	
	mm_prob = sum(mm_freqs.values())

	dr_mm = mm_round(dr_mm_freqs)
	dq_mm = mm_round(dq_mm_freqs)
	dr_mm = str(dr_mm).replace("{","").replace("}", "").replace("defaultdict(<class 'dict'>,", "").replace(")", "")
	dq_mm = str(dq_mm).replace("{","").replace("}", "").replace("defaultdict(<class 'dict'>,", "").replace(")", "")
	#print(mm_freqs)
	dr_mm_freqs = dict( sorted(dr_mm_freqs.items(), key=lambda kv: kv[1], reverse=True))
	dq_mm_freqs = dict( sorted(dq_mm_freqs.items(), key=lambda kv: kv[1], reverse=True))
	mm_freqs = dict( sorted(mm_freqs.items(), key=lambda kv: kv[1], reverse=True))


	mm_prob = round(mm_prob,3)
	return mm_prob,dr_mm, dq_mm
	#return mm_prob,dr_mm_freqs, dq_mm_freqs,ddr1alleles_len,ddr1_cdf_len, ddr2alleles_len,ddr2_cdf_len,ddr_len,ddr_cdf_len,ddq1alleles_len,ddq1_cdf_len,ddq2alleles_len,ddq2_cdf_len,ddq_len,ddq_cdf_len,donor_len,donor_cdf_len, rdr1alleles_len,rdr1_cdf_len,rdr2alleles_len,rdr2_cdf_len,rdr_len,rdr_cdf_len,rdq1alleles_len,rdq1_cdf_len,rdq2alleles_len,rdq2_cdf_len,rdq_len,rdq_cdf_len,recip_len,recip_cdf_len,pairs_len
	#return mm_prob,dr_mm_freqs, dq_mm_freqs,ddr1alleles_len,ddr1_aa_len, ddr1_cdf_len, ddr2alleles_len,ddr2_aa_len,ddr2_cdf_len,ddr_len,ddr_cdf_len,ddq1alleles_len,ddq1_aa_len,ddq1_cdf_len,ddq2alleles_len,ddq2_aa_len,ddq2_cdf_len,ddq_len,ddq_cdf_len,donor_len,donor_cdf_len, rdr1alleles_len,rdr1_aa_len,rdr1_cdf_len,rdr2alleles_len,rdr2_aa_len,rdr2_cdf_len,rdr_len,rdr_cdf_len,rdq1alleles_len,rdq1_aa_len,rdq1_cdf_len,rdq2alleles_len,rdq2_aa_len,rdq2_cdf_len,rdq_len,rdq_cdf_len,recip_len,recip_cdf_len,pairs_len

def DRDQFIBERSE_Alleles(ddr1alleles,ddr2alleles, rdr1alleles, rdr2alleles, ddq1alleles,ddq2alleles,rdq1alleles,rdq2alleles):
	#Dictionaries to store D or R Haplotype Frequencies
	#print("Memory when initiating AA assignment and Match Function:", psutil.virtual_memory())
	#start = timeit.default_timer()
	start_all = timeit.default_timer()
	#Generate dicts to calculate allele freqs to probs and use probs throughout
	ddr1_probs = defaultdict(dict)
	ddr2_probs = defaultdict(dict)
	ddq1_probs = defaultdict(dict)
	ddq2_probs = defaultdict(dict)

	rdr1_probs = defaultdict(dict)
	rdr2_probs = defaultdict(dict)
	rdq1_probs = defaultdict(dict)
	rdq2_probs = defaultdict(dict)

	#individual dicts for single locus and multi loci D and R cals
	dr_list_donor = defaultdict(dict)
	dq_list_donor = defaultdict(dict)
	geno_list_donor = defaultdict(dict)

	dr_list_recip = defaultdict(dict)
	dq_list_recip = defaultdict(dict)
	geno_list_recip = defaultdict(dict)

	#Dictionary to store D|R Haplotype pair frequency
	db_geno = defaultdict(dict)

	#Dictionary to store D|R Haplotype pair probability
	#db_prob= defaultdict(dict)

	#Set AA string length
	#dr_start_position = ard_start_pos[drloc]
	#dr_end_position = ard_end_pos[drloc]
	
	#dq_start_position = ard_start_pos[dqloc]
	#dq_end_position = ard_end_pos[dqloc]

	#Dictionaries to store AA-MM frequencies not MM or MM
	dr_pos_freqs = defaultdict(dict)
	dq_pos_freqs = defaultdict(dict)

	dr_mm_freqs = defaultdict(dict)
	dq_mm_freqs = defaultdict(dict)

	#this will give haplotype D|R pair AA-MM prob
	#all_pos_freqs = defaultdict(dict)
	mm_freqs = defaultdict(dict)
	match_dict = defaultdict(dict)

	#calculate donor and recip DR and DQ frequencies and then create a genotype dictionary for each
	#calculate DR loci freq from allele level freqs
	#sort dictionary ascending and apply CDF function to cut low probabilities for optomization

	total_recip_dr1_freq = sum(rdr1alleles.values())
	total_recip_dr2_freq = sum(rdr2alleles.values())
	total_donor_dr1_freq = sum(ddr1alleles.values())
	total_donor_dr2_freq = sum(ddr2alleles.values())

	total_recip_dq1_freq = sum(rdq1alleles.values())
	total_recip_dq2_freq = sum(rdq2alleles.values())
	total_donor_dq1_freq = sum(ddq1alleles.values())
	total_donor_dq2_freq = sum(ddq2alleles.values())

	for key, value in ddr1alleles.items():
		ddr1_probs[key] = (value/total_donor_dr1_freq)
	for key, value in ddr2alleles.items():
		ddr2_probs[key] = (value/total_donor_dr2_freq)
	for key, value in ddq1alleles.items():
		ddq1_probs[key] = (value/total_donor_dq1_freq)
	for key, value in ddq2alleles.items():
		ddq2_probs[key] = (value/total_donor_dq2_freq)


	for key, value in rdr1alleles.items():
		rdr1_probs[key] = (value/total_recip_dr1_freq)
	for key, value in rdr2alleles.items():
		rdr2_probs[key] = (value/total_recip_dr2_freq)
	for key, value in rdq1alleles.items():
		rdq1_probs[key] = (value/total_recip_dq1_freq)
	for key, value in rdq2alleles.items():
		rdq2_probs[key] = (value/total_recip_dq2_freq)

	#print(ddr1_probs)
	#ddr1_cdf = cumulative(ddr1alleles)
	#ddr2_cdf = cumulative(ddr2alleles)
	ddr1_cdf = cumulative(ddr1_probs)
	ddr2_cdf = cumulative(ddr2_probs)
	for ddra1 in ddr1_cdf:
		for ddra2 in ddr2_cdf:
			#ddra = (ddra1)+ '+' + (ddra2)
			#da_un =ddra1 + ddra2
			#ddra = '+'.join(sorted(da_un, key=str.lower))
			ddra = '+'.join(sort_string([ddra1,ddra2]))
			#print(ddra)
			if (ddra1 == ddra2):
				ddrf1 = ddr1alleles.get(ddra1)
				gdrfreq = ddrf1*ddrf1
				if ddra not in dr_list_donor:
					dr_list_donor[ddra]=gdrfreq
				else:
					dr_list_donor[ddra]=dr_list_donor[ddra] + gdrfreq
			else:
				df1 = ddr1alleles.get(ddra1)
				df2 = ddr2alleles.get(ddra2)
				gdrfreq = 2*df1*df2
				if ddra not in dr_list_donor:
					dr_list_donor[ddra]=gdrfreq
				else:
					dr_list_donor[ddra]=dr_list_donor[ddra] + gdrfreq

	#ddq1_cdf = cumulative(ddq1alleles)
	#ddq2_cdf = cumulative(ddq2alleles)
	ddq1_cdf = cumulative(ddq1_probs)
	ddq2_cdf = cumulative(ddq2_probs)

	#calculate DQ loci freq from allele level freqs
	for ddqa1 in ddq1_cdf:
		for ddqa2 in ddq2_cdf:
			#ddqa = (ddqa1)+ '+' + (ddqa2)
			ddqa = '+'.join(sort_string([ddqa1,ddqa2]))
			#daq_un =ddqa1 + ddqa2
			#ddqa = '+'.join(sorted(daq_un, key=str.lower))
			if (ddqa1 == ddqa2):
				ddqf1 = ddq1alleles.get(ddqa1)
				gdqfreq = ddqf1*ddqf1
				if ddqa not in dq_list_donor:
					dq_list_donor[ddqa]=gdqfreq
				else:
					dq_list_donor[ddqa]=dq_list_donor[ddqa] + gdqfreq
			else:
				df1 = ddq1alleles.get(ddqa1)
				df2 = ddq2alleles.get(ddqa2)
				gdqfreq = 2*df1*df2
				if ddqa not in dq_list_donor:
					dq_list_donor[ddqa]=gdqfreq
				else:
					dq_list_donor[ddqa]=dq_list_donor[ddqa] + gdqfreq
	#calculate DR|DQ geno freq from individual DR and DQ loci freqs level freqs
	'''
	for ddr in dr_list_donor:
		for ddq in dq_list_donor:
			dgeno = (ddr)+ '~' + (ddq)
			df1 = dr_list_donor.get(ddr)
			df2 = dq_list_donor.get(ddq)
			gfreq = df1*df2
			if dgeno not in geno_list_donor:
				geno_list_donor[dgeno]=gfreq
			else:
				geno_list_donor[dgeno]=geno_list_donor[dgeno] + gfreq
	'''

	ddr_cdf = cumulative(dr_list_donor)
	ddq_cdf = cumulative(dq_list_donor)

	#gives all key pernumerations and multiplies values
	for key, value in ddr_cdf.items():
		for key2, value2 in ddq_cdf.items():
			geno_list_donor["{}^{}".format(key, key2)] = value * value2
	#print(geno_list_donor)
	#print("Memory after Donor calculations and pernumeration of multi-loci:", psutil.virtual_memory())


	#print("Seconds to calculate Donor DR|DQ Genotype frequency:", timeit.default_timer()-start)
	#start = timeit.default_timer()
	#rdr1_cdf = cumulative(rdr1alleles)
	#rdr2_cdf = cumulative(rdr2alleles)

	rdr1_cdf = cumulative(rdr1_probs)
	rdr2_cdf = cumulative(rdr2_probs)

	#calculate DR loci freq from allele level freqs
	for rdra1 in rdr1_cdf:
		for rdra2 in rdr2_cdf:
			#rdra = (rdra1)+ '+' + (rdra2)
			rdra = '+'.join(sort_string([rdra1,rdra2]))
			#rdra_un =rdra1 + rdra2
			#rdra = '+'.join(sorted(rdra_un, key=str.lower))
			if (rdra1 == rdra2):
				rdrf1 = rdr1alleles.get(rdra1)
				gdrfreq = rdrf1*rdrf1
				if rdra not in dr_list_recip:
					dr_list_recip[rdra]=gdrfreq
				else:
					dr_list_recip[rdra]=dr_list_recip[rdra] + gdrfreq
			else:
				df1 = rdr1alleles.get(rdra1)
				df2 = rdr2alleles.get(rdra2)
				gdrfreq = 2*df1*df2
				if rdra not in dr_list_recip:
					dr_list_recip[rdra]=gdrfreq
				else:
					dr_list_recip[rdra]=dr_list_recip[rdra] + gdrfreq

	#rdq1_cdf = cumulative(rdq1alleles)
	#rdq2_cdf = cumulative(rdq2alleles)
	rdq1_cdf = cumulative(rdq1_probs)
	rdq2_cdf = cumulative(rdq2_probs)
	#calculate DQ loci freq from allele level freqs
	for rdqa1 in rdq1_cdf:
		for rdqa2 in rdq2_cdf:
			#rdqa = (rdqa1)+ '+' + (rdqa2)
			rdqa = '+'.join(sort_string([rdqa1,rdqa2]))
			#rdaq_un = list(rdqa1) + list(rdqa2)
			#print(rdaq_un)
			#rdqa = '+'.join(sorted(rdaq_un))
			#print(rdqa)
			#print(type(rdqa1))
			#rdqa = '+'.join(sorted(rdqa1,rdqa2 key=str.lower))
			if (rdqa1 == rdqa2):
				rdqf1 = rdq1alleles.get(rdqa1)
				gdqfreq = rdqf1*rdqf1
				if rdqa not in dq_list_recip:
					dq_list_recip[rdqa]=gdqfreq
				else:
					dq_list_recip[rdqa]=dq_list_recip[rdqa] + gdqfreq
			else:
				df1 = rdq1alleles.get(rdqa1)
				df2 = rdq2alleles.get(rdqa2)
				gdqfreq = 2*df1*df2
				if rdqa not in dq_list_recip:
					dq_list_recip[rdqa]=gdqfreq
				else:
					dq_list_recip[rdqa]=dq_list_recip[rdqa] + gdqfreq
					
	#calculate DR|DQ geno freq from individual DR and DQ loci freqs level freqs
	#geno_list_recip= combine_dicts(a, b, operator.mul)	
	#print(geno_list_recip)
	#gives all key pernumerations and multiplies values
	rdr_cdf = cumulative(dr_list_recip)
	rdq_cdf = cumulative(dq_list_recip)
	for key, value in rdr_cdf.items():
		for key2, value2 in rdq_cdf.items():
			geno_list_recip["{}^{}".format(key, key2)] = value * value2
	#print(geno_list_recip)

	'''
	for rdr in dr_list_recip:
		for rdq in dq_list_recip:
			rgeno = (rdr)+ '~' + (rdq)
			df1 = dr_list_recip.get(rdr)
			df2 = dq_list_recip.get(rdq)
			gfreq = df1*df2
			if rgeno not in geno_list_recip:
				geno_list_recip[rgeno]=gfreq
			else:
				geno_list_recip[rgeno]=geno_list_recip[rgeno] + gfreq	
	'''
	#print("Memory after Recipient calculations and pernumeration of multi-loci:", psutil.virtual_memory())
	#print("Seconds to calculate Recipient DR|DQ Genotype frequencies:", timeit.default_timer()-start)
	#start = timeit.default_timer()
	donor_cdf = cumulative(geno_list_donor)
	recip_cdf = cumulative(geno_list_recip)
	####Adding to base off genotype probs instead of using frequencies as requested by Dr.Gragert
	total_donor_geno_freq = sum(donor_cdf.values())
	total_recip_geno_freq = sum(recip_cdf.values())

	#Get strings for donor and recip geno types, calculate pair frequency, and store for probability calculation
	#gives all key pernumerations and multiplies probs
	for key, value in donor_cdf.items():
		for key2, value2 in recip_cdf.items():
			db_geno["{}|{}".format(key, key2)] = (value/total_recip_geno_freq) * (value2/total_donor_geno_freq)
			#db_geno["{}|{}".format(key, key2)] = (value) * (value2)

										
	#print("Memory after D|R pair probability calculations and pernumerations:", psutil.virtual_memory())
	#print("Seconds to calculate D|R Pair Probabilities and Pernumerate Dicts:", timeit.default_timer()-start)
	#start = timeit.default_timer()
	'''
	for dgeno in geno_list_donor:
		for rgeno in geno_list_recip:
			combo = (dgeno) + '|' + (rgeno)
			rfgeno = geno_list_recip.get(rgeno)
			dfgeno = geno_list_donor.get(dgeno)

			rpgeno = rfgeno/total_recip_geno_freq
			dpgeno = dfgeno/total_donor_geno_freq

			# Calculate Pair Probability
			drprob = rpgeno*dpgeno
			if combo not in db_geno:
				db_geno[combo] = drprob
			#else:
				#db_geno[combo] = db_geno[combo] + drprob

			#print("Time to get D|R Pair Probability:", timeit.default_timer()-start)
			#start = timeit.default_timer()
			#Split and prep haplos for AA pos probability by summing D|R pair probability for each possible AA combination of all 4 alleles
			#First get AA assignments for strings based on alleles
	'''
	#db_cdf = cumulative(db_geno)

	for combo in db_geno:
		drprob = db_geno.get(combo)
		dgeno,rgeno = combo.split('|')
		#print(dgeno)
		ddr,ddq = dgeno.split('^')
		#print(ddr)
		da1,da2 = ddr.split('+')
		dq1,dq2 = ddq.split('+')

		dastring1 = drb_string_dict.get(da1)#aa_mm.HLA_seq[da1].seq[dr_start_position-1:dr_end_position]
		dastring2 = drb_string_dict.get(da2)#aa_mm.HLA_seq[da2].seq[dr_start_position-1:dr_end_position]
		dastring1 = str(dastring1)
		dastring2 = str(dastring2)

		dqstring1 = dqb_string_dict.get(dq1)#aa_mm.HLA_seq[dq1].seq[dq_start_position-1:dq_end_position]
		dqstring2 = dqb_string_dict.get(dq2)#aa_mm.HLA_seq[dq2].seq[dq_start_position-1:dq_end_position]
		dqstring1 = str(dqstring1)
		dqstring2 = str(dqstring2)

		rdr,rdq = rgeno.split('^')

		ra1,ra2 = rdr.split('+')
		rq1,rq2 = rdq.split('+')


		rastring1 = drb_string_dict.get(ra1)#aa_mm.HLA_seq[ra1].seq[dr_start_position-1:dr_end_position]
		rastring2 = drb_string_dict.get(ra2)#aa_mm.HLA_seq[ra2].seq[dr_start_position-1:dr_end_position]
		rastring1 = str(rastring1)
		rastring2 = str(rastring2)

		rqstring1 = dqb_string_dict.get(rq1)#aa_mm.HLA_seq[rq1].seq[dq_start_position-1:dq_end_position]
		rqstring2 = dqb_string_dict.get(rq2)#aa_mm.HLA_seq[rq2].seq[dq_start_position-1:dq_end_position]
		rqstring1 = str(rqstring1)
		rqstring2 = str(rqstring2)
		#print("Time to get DR and DQ AA strings:", timeit.default_timer()-start)

		#start = timeit.default_timer()
		#Next based on loci, DR or DQ, check for MM at select positions set at the begining of the script
		for pos in list(dr_fibers_pos):
			daa1 = dastring1[pos-1]
			daa2 = dastring2[pos-1]
			raa1 = rastring1[pos-1]
			raa2 = rastring2[pos-1]
			
			#sort the AA asignmnets of D and R alphabetically to compare
			daa = daa1 + daa2
			daa = ''.join(sorted(daa, key=str.lower))
			raa = raa1 + raa2
			raa = ''.join(sorted(raa, key=str.lower))

			# join sorted AA assignments in D|R formation for comparison
			show = daa + '|'+ raa

			#Check for homo AA asignment in Donor
			if (daa1 == daa2):
				#if Donor is Homo and Matches at least 1/2 Recip AA it is a match due to TX directionality
				if (daa1 == raa1 or daa1==raa2):
					#continue
					#store D|R AA assignments and probability of Matches
					if pos not in match_dict["Match"]:
						match_dict["Match"][pos]=drprob
					else:
						match_dict["Match"][pos]=match_dict["Match"][pos]+drprob
					if show not in dr_pos_freqs[pos]:
						dr_pos_freqs[pos][show] = drprob
					else: 
						dr_pos_freqs[pos][show] = dr_pos_freqs[pos][show] + drprob
				else:
					#store D|R AA assignments and probability of MisMatches
					if show not in dr_mm_freqs[pos]:
						dr_mm_freqs[pos][show] = drprob
					else: 
						dr_mm_freqs[pos][show] = dr_mm_freqs[pos][show] + drprob

					if pos not in match_dict["MM"]:
						match_dict["MM"][pos]=drprob
					else:
						match_dict["MM"][pos]=match_dict["MM"][pos]+drprob
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue
			#Else statement if Donor AA assignment is heterozygous
			else:		
				if (raa == daa):
					#continue
					#store D|R AA assignments and probability of Matches
					if show not in dr_pos_freqs[pos]:
						dr_pos_freqs[pos][show] = drprob
					else: 
						dr_pos_freqs[pos][show] = dr_pos_freqs[pos][show] + drprob
					if pos not in match_dict["Match"]:
						match_dict["Match"][pos]=drprob
					else:
						match_dict["Match"][pos]=match_dict["Match"][pos]+drprob

				else:
					#store D|R AA assignments and probability of Mismatches
					if show not in dr_mm_freqs[pos]:
						dr_mm_freqs[pos][show] = drprob
					else: 
						dr_mm_freqs[pos][show] = dr_mm_freqs[pos][show] + drprob
					if pos not in match_dict["MM"]:
						match_dict["MM"][pos]=drprob
					else:
						match_dict["MM"][pos]=match_dict["MM"][pos]+drprob					
					#store D|R genotype assignments and probability of those that MM to get MM probability for DR and DQ overall
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue

		#print("Time to get DR AA MM positonal probabilities:", timeit.default_timer()-start)
		#start = timeit.default_timer()

		#DR methodology repeated for DQ positons 
		for pos in list(dq_fibers_pos):
			dqaa1 = dqstring1[pos-1]
			dqaa2 = dqstring2[pos-1]
			rqaa1 = rqstring1[pos-1]
			rqaa2 = rqstring2[pos-1]
			dqaa = dqaa1 + dqaa2
			dqaa = ''.join(sorted(dqaa, key=str.lower))
			rqaa = rqaa1 + rqaa2
			rqaa = ''.join(sorted(rqaa, key=str.lower))
			show = dqaa + '|'+ rqaa
			if (dqaa1==dqaa2):
				if (dqaa1==rqaa1 or dqaa1==rqaa2):
					#continue
					if show not in dq_pos_freqs[pos]:
						dq_pos_freqs[pos][show] = drprob
					else: 
						dq_pos_freqs[pos][show] = dq_pos_freqs[pos][show] + drprob
					if pos not in match_dict["Match"]:
						match_dict["Match"][pos]=drprob
					else:
						match_dict["Match"][pos]=match_dict["Match"][pos]+drprob
				else:
					if show not in dq_mm_freqs[pos]:
						dq_mm_freqs[pos][show] = drprob
					else: 
						dq_mm_freqs[pos][show] = dq_mm_freqs[pos][show] + drprob
					if pos not in match_dict["MM"]:
						match_dict["MM"][pos]=drprob
					else:
						match_dict["MM"][pos]=match_dict["MM"][pos]+drprob
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue  
			else:
				if (rqaa == dqaa):
					#continue
					if show not in dq_pos_freqs[pos]:
						dq_pos_freqs[pos][show] = drprob
					else:
						dq_pos_freqs[pos][show] = dq_pos_freqs[pos][show] + drprob
					if pos not in match_dict["Match"]:
						match_dict["Match"][pos]=drprob
					else:
						match_dict["Match"][pos]=match_dict["Match"][pos]+drprob
				else:
					if show not in dq_mm_freqs[pos]:
						dq_mm_freqs[pos][show] = drprob
					else: 
						dq_mm_freqs[pos][show] = dq_mm_freqs[pos][show] + drprob
					if pos not in match_dict["MM"]:
						match_dict["MM"][pos]=drprob
					else:
						match_dict["MM"][pos]=match_dict["MM"][pos]+drprob
					if combo not in mm_freqs:
						mm_freqs[combo] = drprob
					else:
						continue 
			#else:
				#continue	
	#print("Memory after D|R pair AA probability summation:", psutil.virtual_memory())
	#print("Seconds to get DR and DQ AA probabilities:", timeit.default_timer()-start)

	#mm_freqs = dict( sorted(mm_freqs.items(), key=lambda kv: kv[1], reverse=False))
	#db_geno = dict( sorted(db_geno.items(), key=lambda kv: kv[1], reverse=False))



	#ddq1alleles = dict( sorted(ddq1alleles.items(), key=lambda kv: kv[1], reverse=False))
	#ddq2alleles = dict( sorted(ddq2alleles.items(), key=lambda kv: kv[1], reverse=False))
	#dr_list_donor = dict( sorted(dr_list_donor.items(), key=lambda kv: kv[1], reverse=False))
	#dq_list_donor = dict( sorted(dq_list_donor.items(), key=lambda kv: kv[1], reverse=False))
	#geno_list_donor = dict( sorted(geno_list_donor.items(), key=lambda kv: kv[1], reverse=False))

	#rdr1alleles = dict( sorted(rdr1alleles.items(), key=lambda kv: kv[1], reverse=False))
	#rdr2alleles = dict( sorted(rdr2alleles.items(), key=lambda kv: kv[1], reverse=False))
	#rdq1alleles = dict( sorted(rdq1alleles.items(), key=lambda kv: kv[1], reverse=False))
	#rdq2alleles = dict( sorted(rdq2alleles.items(), key=lambda kv: kv[1], reverse=False))
	#dr_list_recip = dict( sorted(dr_list_recip.items(), key=lambda kv: kv[1], reverse=False))
	#dq_list_recip = dict( sorted(dq_list_recip.items(), key=lambda kv: kv[1], reverse=False))
	#geno_list_recip = dict( sorted(geno_list_recip.items(), key=lambda kv: kv[1], reverse=False))

	#vals = np.fromiter(ddr1alleles.values(), dtype=float)
	#vals_sorted = np.sort(vals)
	#count, bins_count = np.histogram(vals, bins=len(vals))
	#pdf = count / sum(count)
  
	# using numpy np.cumsum to calculate the CDF
	# We can also find using the PDF values by looping and adding
	#cdf = np.cumsum(pdf)
	
	# plotting PDF and CDF
	#plt.plot(bins_count[1:], pdf, color="red", label="PDF")
	#plt.show()
	#plt.plot(bins_count[1:], cdf, label="CDF")
	#plt.legend()
	#plt.show()


	# calculate the proportional values of samples
	#https://stackoverflow.com/questions/24788200/calculate-the-cumulative-distribution-function-cdf-in-python
	#doesnt work gives linear values
	#p = 1. * np.arange(len(vals)) / (len(vals) - 1)
	#fig = plt.figure()
	##ax2 = fig.add_subplot(122)
	##ax2.plot(p, vals_sorted)
	#ax2.set_xlabel('$x$')
	#ax2.set_ylabel('$p$')
	#plt.show()

	#print("Memory after reordering dictionaires:", psutil.virtual_memory())

	#print("Dictionary of Donor DR Antigen 1 Alleles and Frequencies:",ddr1alleles)
	#ddr1 = cumulative(ddr1alleles)
	#print("Dictionary of Donor DR Antigen 1 Alleles and Frequencies Cumulative Distribution:",ddr1)
	#print("Dictionary of Donor DR Antigen 2 Alleles and Frequencies:",ddr2alleles)
	#ddr2 = cumulative(ddr2alleles)
	#print("Dictionary of Donor DR Antigen 2 Alleles and Frequencies Cumulative Distribution:",ddr2)
	#print("Dictionary of Donor DR loci frequencies calculated from allele level frequencies:",dr_list_donor)
	#ddr = cumulative(dr_list_donor)
	#print("Dictionary of Donor DR loci frequencies calculated from allele level frequencies Cumulative Distribution:",ddr)
	#print("Dictionary of Donor DQ Antigen 1 Alleles and Frequencies:",ddq1alleles)
	#ddq1 = cumulative(ddq1alleles)
	#print("Dictionary of Donor DQ Antigen 1 Alleles and Frequencies Cumulative Distribution:",ddq1)
	#print("Dictionary of Donor DQ Antigen 2 Alleles and Frequencies:",ddq2alleles)
	#ddq2 = cumulative(ddq2alleles)
	#print("Dictionary of Donor DQ Antigen 2 Alleles and Frequencies Cumulative Distribution:",ddq2)
	#print("Dictionary of Donor DQ loci frequencies calculated from allele level frequencies:", dq_list_donor)
	#dq =cumulative(dq_list_donor)
	#print("Dictionary of Donor DQ loci frequencies calculated from allele level frequencies Cumulative Distribution:", dq)	
	#print("Dictionary of Donor DR_DQ Genotype frequencies calculated from loci level frequencies:", geno_list_donor)
	#geno = cumulative(geno_list_donor)
	#print("Dictionary of Donor DR_DQ Genotype frequencies calculated from loci level frequencies Cumulative Distribution:", geno)

	#print("Dictionary of Recip DR Antigen 1 Alleles and Frequencies:",rdr1alleles)
	#ddr1 = cumulative(rdr1alleles)
	#print("Dictionary of Recip DR Antigen 1 Alleles and Frequencies Cumulative Distribution:",ddr1)
	#print("Dictionary of Recip DR Antigen 2 Alleles and Frequencies:",rdr2alleles)
	#ddr2 = cumulative(rdr2alleles)
	#print("Dictionary of Recip DR Antigen 2 Alleles and Frequencies Cumulative Distribution:",ddr2)
	#print("Dictionary of Recip DR loci frequencies calculated from allele level frequencies:",dr_list_recip)
	#ddr = cumulative(dr_list_recip)
	#print("Dictionary of Recip DR loci frequencies calculated from allele level frequencies Cumulative Distribution:",ddr)
	#print("Dictionary of Recip DQ Antigen 1 Alleles and Frequencies:",rdq1alleles)
	#ddq1 = cumulative(rdq1alleles)
	#print("Dictionary of Recip DQ Antigen 1 Alleles and Frequencies Cumulative Distribution:",ddq1)
	#print("Dictionary of Recip DQ Antigen 2 Alleles and Frequencies:",rdq2alleles)
	#ddq2 = cumulative(rdq2alleles)
	#print("Dictionary of Recip DQ Antigen 2 Alleles and Frequencies Cumulative Distribution:",ddq2)	
	#print("Dictionary of Recip DQ loci frequencies calculated from allele level frequencies:", dq_list_recip)
	#dq =cumulative(dq_list_recip)
	#print("Dictionary of Recip DQ loci frequencies calculated from allele level frequencies Cumulative Distribution:", dq)	
	#print("Dictionary of Recip DR_DQ Genotype frequencies calculated from loci level frequencies:", geno_list_recip)
	#geno = cumulative(geno_list_donor)
	#print("Dictionary of Recip DR_DQ Genotype frequencies calculated from loci level frequencies Cumulative Distribution:", geno)
	#print("Loci analyzed for High Risk AA-MM Pos:", locus)
	#print("Dictionary of Donor DR loci frequencies calculated from allele level frequencies:", dict(islice(dr_list_donor.items(), 0, 4)))
	#print("Dictionary of Donor DQ loci frequencies calculated from allele level frequencies:", dict(islice(dq_list_donor.items(), 0, 4)))	
	#print("Dictionary of Donor DR_DQ Genotype frequencies calculated from loci level frequencies:", dict(islice(geno_list_donor.items(), 0, 4)))
	#print("Total number of Donor DR_DQ types in Dictionary:",len(geno_list_donor))
	#print("Dictionary of Recip DR_DQ Genotype frequencies calculated from loci level frequencies:", dict(islice(geno_list_recip.items(), 0, 4)))
	#print("Total number of Recip DR_DQ types in Dictionary:",len(geno_list_recip))

	#print("Dictionary of D|R  DR_DQ Genotype pair frequencies calculated from DR_DQ D|R Genotype frequencies:", dict(islice(db_geno.items(), 0, 100)))
	#print("Dictionary of D|R  DR_DQ Genotype pair frequencies calculated from DR_DQ D|R Genotype frequencies:", db_geno)
	#geno_prob = round(sum(db_geno.values()),2)
	geno_prob = sum(db_geno.values())
	#print("Sum of Pair Probability Dictionary:",geno_prob)
	#print("Total number of D|R Pairs in Dictionary:",len(db_geno))	

	#print("Dictionary of all DR AA freqs for D|R haplotype pairs:", dict(islice(dr_pos_freqs.items(), 0, 4)))
	#print("Dictionary of all DR AA-MM freqs of D|R haplotype pairs identified as HighRisk:", dict(islice(dr_mm_freqs.items(), 0, 4)))
	#print("Dictionary of all DR AA Match probs for D|R haplotype pairs:", dict(dr_pos_freqs))
	#print("Dictionary of all DQ AA Match probs of D|R haplotype pairs :", dict(dq_pos_freqs))
	#print("Dictionary of all DR AA-MM probs for D|R haplotype pairs:", dict(dr_mm_freqs))
	#print("Dictionary of all DQ AA-MM probs of D|R haplotype pairs :", dict(dq_mm_freqs))
	
	mm_prob = sum(mm_freqs.values())
	#all_pos = sum(all_pos_freqs.values())
	#mm_prob = round(sum(mm_freqs.values()),2)

	#db_prob= sum(db_geno.values())

	#print("Dictionary of D|R  DR_DQ Genotype pair probabilities calculated from DR_DQ D|R Genotype frequencies:", dict(islice(db_geno.items(), 0, 4)))
	#print("Total number of D|R Pairs in Dictionary:",len(db_geno))	
	#print("Total of probability in Dictionary:",db_prob)

	#print("Dictionary of D|R  DR_DQ Genotype pair probabilities calculated from DR_DQ D|R Genotype frequencies with selected positions that MM:", dict(islice(mm_freqs.items(), 0, 100)))
	#print("Dictionary of D|R  DR_DQ Genotype pair probabilities calculated from DR_DQ D|R Genotype frequencies with selected positions that MM:", mm_freqs)
	#print("Sum of Pair Probability in MM Dictionary:",mm_prob)
	#print("Total number of D|R Pairs in MM Dictionary:",len(mm_freqs))

	#print("Dictionary of D|R  DR_DQ Genotype pair probabilities calculated from DR_DQ D|R Genotype frequencies with >=1 AA-MM:", dict(islice(mm_freqs.items(), 0, 4)))
	#print("Probability of D|R pair having >= 1 AA-MM Identified as FIBERS High Risk:", mm_prob)
	#dr_sum = sum(dr_mm_freqs.values())
	#dq_sum = sum(dq_mm_freqs.values())
	#pos_sum = dr_sum + dq_sum
	#print(pos_sum)
	#print("Probability of D|R pair having >= 1 AA-MM Identified as FIBERS High Risk:", pos_sum)

	#print("Dictionary of all DR AA Match probs for D|R haplotype pairs:", dict(dr_pos_freqs))
	#print("Dictionary of all DQ AA Match probs of D|R haplotype pairs :", dict(dq_pos_freqs))

	#print("Dictionary of all DR AA-MM probs for D|R haplotype pairs:", dict(dr_mm_freqs))
	#print("Dictionary of all DQ AA-MM probs of D|R haplotype pairs :", dict(dq_mm_freqs))
	
	#print("Dictionary seperated by Match vs MM:", dict(match_dict))

	res = defaultdict(int)
	for inner_dict in match_dict.values():
		for key, value in inner_dict.items():
			res[key] += value
 
	# converting defaultdict to dictionary
	res = dict(res)
	#print("Summation of Mismatch vs Match dictionary by AA position is : " + str(res))
	#print("Seconds for Duration of Entire Function:", timeit.default_timer()-start_all)
	#dr_mm=format_mm(dr_mm_freqs)
	#dq_mm=format_mm(dq_mm_freqs)
	dr_mm = mm_round(dr_mm_freqs)
	dq_mm = mm_round(dq_mm_freqs)
	dr_mm = str(dr_mm).replace("{","").replace("}", "").replace("defaultdict(<class 'dict'>,", "").replace(")", "").replace("'", "")
	dq_mm = str(dq_mm).replace("{","").replace("}", "").replace("defaultdict(<class 'dict'>,", "").replace(")", "").replace("'", "")
	#dr_mm = round(dr_mm,2)
	return geno_prob,mm_prob,dr_mm, dq_mm

# weighted choice from https://scaron.info/blog/python-weighted-choice.html
def weighted_choice(seq, weights):
    assert len(weights) == len(seq)
    assert abs(1. - sum(weights)) < 1e-6

    x = random.random()
    for i, elmt in enumerate(seq):
        if x <= weights[i]:
            return elmt
        x -= weights[i]
