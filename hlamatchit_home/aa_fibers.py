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
import pandas as pd
import random
from .aa_matching_msf import *
aa_mm = AAMatch(dbversion=3500)
import re




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



def afibershazard(pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.split(",")
	ahazard = 0
	pos2_list = []
	for pos in pos_list:
		if (pos == "12") or (pos =="44") or (pos=="63") or (pos=="105") or (pos == "111") or (pos=="114") or (pos=="152") or (pos=="161") or (pos == "166") or (pos== "167"):
			ahazard = 1.09
			pos2_list.append(pos)
	pos2_list = ', '.join(str(item) for item in pos2_list)
	return ahazard, pos2_list




def cfibershazard(pos_list):
#print(locus)
#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.split(",")
	chazard = 0
	pos2_list = []
	for pos in pos_list:
		if (pos== "11") or (pos== "35") or (pos== "30") or (pos=="31") or (pos=="39") or (pos=="41") or (pos=="46") or (pos== "65") or (pos=="70") or (pos=="97") or (pos=="108") or (pos=="122") or (pos=="143") or (pos=="156") or (pos=="160") or (pos=="163") or (pos=="176") or (pos=="179"):
			chazard = 1.04
			pos2_list.append(pos)
	pos2_list = ', '.join(str(item) for item in pos2_list)
	return chazard, pos2_list




def bfibershazard(pos_list):
#print(locus)
#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.split(",")
	bhazard = 0
	pos2_list = []
	for pos in pos_list:
		if (pos == "23") or (pos=="24") or (pos=="46") or (pos=="67") or (pos=="136") or (pos=="145"):
			bhazard = 1.04
			pos2_list.append(pos)
	pos2_list = ', '.join(str(item) for item in pos2_list)
	return bhazard, pos2_list



def drfibershazard(pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.split(",")
	drhazard = 0
	pos2_list = []
	prob = 0
	for pos in pos_list:
		if (pos== "11") or (pos == "14") or (pos=="16") or (pos=="23") or (pos=="26") or (pos=="28") or (pos=="30") or (pos=="32") or (pos=="37") or (pos=="50") or (pos=="51") or (pos=="60") or (pos== "78"):
			drhazard = 1.11
			pos2_list.append(pos)
		if (pos == "13") or (pos == "26"):
			prob = 1
		else:
			continue
	return drhazard, pos2_list, prob


def dqfibershazard(pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.split(",")
	dqhazard = 0
	pos2_list = []
	prob = 0
	for pos in pos_list:
		if (pos == "18") or (pos == "49") or (pos == "55") or (pos == "66") or (pos == "74"):
			dqhazard = 1.07
			pos2_list.append(pos)
		if (pos == "30") or (pos == "55"):
			prob = 1
		else:
			continue
	return dqhazard, pos2_list, prob


def drfibersprob(pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.replace("'","")
	pos_list = pos_list.split(",")
	pos2_list = []
	for pos in pos_list:
		#print(pos,type(pos))
		if (pos == "9") or (pos == "10") or (pos == "11") or (pos== "12") or (pos == "13") or (pos == "26") or (pos == "28") or (pos == "30"):
			pos2_list.append(pos)
			print("AA-MM:", pos)
		else:
			continue
	pos2_list = ', '.join(str(item) for item in pos2_list)
	cpos = len(pos2_list)
	if (cpos > 0):
		prob = 1
	else:
		prob = 0
	return pos2_list, prob


def dqfibersprob(pos_list):
	#print(locus)
	#print(pos_list)
	pos_list = pos_list.replace(" ","")
	pos_list = pos_list.replace("'","")
	pos_list = pos_list.split(",")
	pos2_list = []
	for pos in pos_list:
		#print(pos,type(pos))
		if (pos == "30") or (pos == "38") or (pos=="53") or (pos == "55") or (pos == "66") or (pos== "67") or (pos== "71") or (pos == "74") or (pos== "77") or (pos== "84") or (pos == "85") or (pos=="89") or (pos== "90"):
			pos2_list.append(pos)
		else:
			continue
	pos2_list = ', '.join(str(item) for item in pos2_list)
	cpos = len(pos2_list)
	if (cpos > 0):
		prob = 1
	else:
		prob = 0
	return pos2_list, prob


def getAAgenostringmatchFIBERSDQ(dalleles1, dalleles2, ralleles1, ralleles2, locus):
	probs = defaultdict(dict)
	geno_list = defaultdict(dict)
	for da1 in dalleles1:
		for da2 in dalleles2:
			for ra1 in ralleles1:
				for ra2 in ralleles2:
					da = (da1)+ '+' + (da2)
					#print(daa)
					#da = '+'.join(sorted(da, key=str.lower))
					#print(daa)
					ra = (ra1) + '+' + (ra2)
					#a = '+'.join(sorted(ra, key=str.lower))
					combo = (da) + '|' + (ra)
					#print("combo:",combo)
					df1 = dalleles1.get(da1)
					df2 = dalleles2.get(da2)
					#ddr1 = df1/d1sumant
					#ddr2 = df2/d2sumant
					#ddr = ddr1 + ddr2
					#ddrprob = ddr/2

					rf1 = ralleles1.get(ra1)
					rf2 = ralleles2.get(ra2)
					#rdr1 = rf1/r1sumant
					#rdr2 = rf2/r2sumant

					#rdr = rdr1 + rdr2
					#rdrprob = rdr/2

					#drprob = ddrprob + rdrprob
					#dfprob = drprob/2
					dfprob = df1 +df2 +rf1 +rf2
					geno_list[combo] = dfprob
					start_position = ard_start_pos[locus]
					end_position = ard_end_pos[locus]
					dstring1 = aa_mm.HLA_seq[da1].seq[start_position-1:end_position]
					dstring2 = aa_mm.HLA_seq[da2].seq[start_position-1:end_position]
					rstring1 = aa_mm.HLA_seq[ra1].seq[start_position-1:end_position]
					rstring2 = aa_mm.HLA_seq[ra2].seq[start_position-1:end_position]
					dstring1 = str(dstring1)
					dstring2 = str(dstring2)
					rstring1 = str(rstring1)
					rstring2 = str(rstring2)
					start_position = int(start_position)
					end_position = int(end_position)
					for pos in range(start_position,end_position):
						#print(pos)
						if (pos == 30) or (pos == 38) or (pos==53) or (pos == 55) or (pos == 66) or (pos== 67) or (pos== 71) or (pos == 74) or (pos== 77) or (pos== 84) or (pos == 85) or (pos==89) or (pos== 90):
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
							if(raa == daa):
								continue
							else:
								if combo not in probs:
									probs[combo]=dfprob
								else:
									continue
						else:
							continue

	sumgeno = sum(geno_list.values())
	sumpossible =sum(probs.values())
	#nprobs = len(probs)
	#print(nprobs)
	#sumant= rdsumant*nprobs
	#print(sumant)
	#print(sumpossible)
	#mm_prob = sumpossible/nprobs
	mm_prob = sumpossible/sumgeno
	print(mm_prob)

	return mm_prob


def getAAgenostringmatchFIBERSDR(dalleles1, dalleles2, ralleles1, ralleles2, locus):
	probs = defaultdict(dict)
	geno_list = defaultdict(dict)
	for da1 in dalleles1:
		for da2 in dalleles2:
			for ra1 in ralleles1:
				for ra2 in ralleles2:
					da = (da1)+ '+' + (da2)
					#print(daa)
					#da = '+'.join(sorted(da, key=str.lower))
					#print(daa)
					ra = (ra1) + '+' + (ra2)
					#a = '+'.join(sorted(ra, key=str.lower))
					combo = (da) + '|' + (ra)
					#print("combo:",combo)
					df1 = dalleles1.get(da1)
					df2 = dalleles2.get(da2)
					#ddr1 = df1/d1sumant
					#ddr2 = df2/d2sumant
					#ddr = ddr1 + ddr2
					#ddrprob = ddr/2

					rf1 = ralleles1.get(ra1)
					rf2 = ralleles2.get(ra2)
					#rdr1 = rf1/r1sumant
					#rdr2 = rf2/r2sumant

					#rdr = rdr1 + rdr2
					#rdrprob = rdr/2

					#drprob = ddrprob + rdrprob
					#dfprob = drprob/2
					dfprob = df1 +df2 +rf1 +rf2
					geno_list[combo] = dfprob
					start_position = ard_start_pos[locus]
					end_position = ard_end_pos[locus]
					dstring1 = aa_mm.HLA_seq[da1].seq[start_position-1:end_position]
					dstring2 = aa_mm.HLA_seq[da2].seq[start_position-1:end_position]
					rstring1 = aa_mm.HLA_seq[ra1].seq[start_position-1:end_position]
					rstring2 = aa_mm.HLA_seq[ra2].seq[start_position-1:end_position]
					dstring1 = str(dstring1)
					dstring2 = str(dstring2)
					rstring1 = str(rstring1)
					rstring2 = str(rstring2)
					start_position = int(start_position)
					end_position = int(end_position)
					for pos in range(start_position,end_position):
						#print(pos)
						if (pos == 9) or (pos == 10) or (pos == 11) or (pos== 12) or (pos == 13) or (pos == 26) or (pos == 28) or (pos == 30):
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
								print(raa,daa)
								if combo not in probs:
									probs[combo]=dfprob
								else:
									continue
						else:
							continue
	#print(probs)
	sumgeno = sum(geno_list.values())
	sumpossible =sum(probs.values())
	#nprobs = len(probs)
	#print(nprobs)
	#sumant= rdsumant*nprobs
	print(sumgeno)
	print(sumpossible)
	#mm_prob = sumpossible/nprobs
	mm_prob = sumpossible/sumgeno
	print(mm_prob)

	return mm_prob


#functions for enumarating probs as shown in DRDQFIBERS.ipynb

def antigen2HFalleleE(race,antigen):
	alleles = antigen2allele(antigen)
	print( "Race:", race, "antigen:", antigen)
	print("Alleles from OPTN selected by antigen2allele function:", alleles)
	#print(alleles)
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
	print("Only Alleles with frequency in Genetic Reference Population:", alleles)
	print("Sum of Allele Frequencies for Probability Calcs:", sumant)
	print( "Most probable Allele and its Frequency:", max_allele, max_freq)
	allelefreqs = probcheck(alleles, sumant)
	print("Prob check function calculating Allele level Probs using AlleleFreq/AntigenFreq for every possible Allele:", allelefreqs)
	prob = max_freq/sumant
	print("Most Probable Alleles Probability:", prob)
	return max_allele, max_freq, nalleles, sumant, allelefreqs, prob


def probcheck(possible_alleles,sumant):
	probs= []
	for allele in possible_alleles:
		probs.append(allele)
		freq = possible_alleles.get(allele)
		prob = freq/sumant
		probs.append(prob)
	return probs


def getAAgenostringmatchE(da1, da2, ra1, ra2, locus, dprob, rprob):
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
		if (pos == 9) or (pos == 10) or (pos == 11) or (pos== 12) or (pos == 13) or (pos == 26) or (pos == 28) or (pos == 30):
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
				print("getAAgenostringmatchE AA-MM:", pos, raa,"+",daa)
				prob = 1-((rprob + dprob)/2)
				print("Probability,1-(rprob+dprob/2), of AA-MM given matching assignments for the most probable genotypes:", prob)
			else:
				print("getAAgenostringmatchE AA-MM:", pos, raa,"+",daa)
				prob = (rprob + dprob)/2
				print("Probability, rprob+dprob/2, of AA-MM given mismatching assignments for the most probable genotypes:", prob)
				mm_count+=1
				pos_list.append(pos)
				pos1_list.append(pos)
				combo = daa + "|" + raa
				pos1_list.append(combo)
		else:
			continue

	dstring1 = str(dstring1)
	dstring2 = str(dstring2)
	rstring1 = str(rstring1)
	rstring2 = str(rstring2)
	pos_list = ', '.join(str(item) for item in pos_list)
	pos1_list = ', '.join(str(item) for item in pos1_list)
	start_position = int(start_position)
	end_position = int(end_position)
	return dstring1, dstring2, rstring1, rstring2, mm_count, pos_list, pos1_list, end_position

# weighted choice from https://scaron.info/blog/python-weighted-choice.html
def weighted_choice(seq, weights):
    assert len(weights) == len(seq)
    assert abs(1. - sum(weights)) < 1e-6

    x = random.random()
    for i, elmt in enumerate(seq):
        if x <= weights[i]:
            return elmt
        x -= weights[i]
