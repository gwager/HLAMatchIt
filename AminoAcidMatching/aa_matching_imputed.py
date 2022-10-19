#!/usr/bin/env python3

# aa_matching.py - Module for amino acid matching functions

import re
from os import path
from pathlib import Path
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import random
import requests


# fix for filenotfounderror
pathloc = Path(__file__).parent

# wrapping everything in a class to facilitate import.

class AAMatch:

    def __init__(self, dbversion=3450, ungap=True, imputed=False):
        # if ungap = True, use complete reference sequence
        # if ungap = False, use gapped reference sequence
        self.dbversion = dbversion
        self.ungap = ungap
        self.imputed = imputed
        if self.imputed:
            self.fdir = 'ifasta'
            self.ftype = 'fasta'
        else:
            self.fdir = 'msf'
            self.ftype = 'msf'
        hlaProteinOffset = {
            "A" : 0,
            "B" : 0,
            "C" : 0,
            "DRA" : 0,
            "DRB1" : 0,
            "DRB3" : 0,
            "DRB4" : 0,
            "DRB5" : 0,
            "DQA1" : 0,
            "DQB1" : 0,
            "DPA1" : 0,
            "DPB1" : 0,
        }

        first_ten = {
            "A*01:01" : "GSHSMRYFFT",
            "B*07:02" : "GSHSMRYFYT",
            "C*01:02" : "CSHSMKYFFT",
            "DRB1*01:01" : "GDTRPRFLWQ",
            "DRB3*01:01" : "GDTRPRFLEL",
            "DRB4*01:01" : "GDTQPRFLEQ",
            "DRB5*01:01" : "GDTRPRFLQQ",
            "DQA1*01:01" : "EDIVADHVAS",
            "DQB1*05:01" : "RDSPEDFVYQ",
            "DPA1*01:03" : "IKADHVSTYA",
            "DPB1*01:01" : "RATPENYVYQ"
        }

        first_ten_full = {
            "A*01:01:01:01" : "GSHSMRYFFT",
            "B*07:02:01:01" : "GSHSMRYFYT",
            "C*01:02:01:01" : "CSHSMKYFFT",
            "DRB1*01:01:01:01" : "GDTRPRFLWQ",
            "DRB3*01:01:02:01" : "GDTRPRFLEL",
            "DRB4*01:01:01:01" : "GDTQPRFLEQ",
            "DRB5*01:01:01:01" : "GDTRPRFLQQ",
            "DQA1*01:01:01:01" : "EDIVADHVAS",
            "DQB1*05:01:01:01" : "RDSPEDFVYQ",
            "DPA1*01:03:01:01" : "IKADHVSTYA",
            "DPB1*01:01:01:01" : "RATPENYVYQ"
        }

        ard_first_ten = {
            "A*01:01" : "SHSMRYFFTS",
            "B*07:02" : "SHSMRYFYTS",
            "C*01:02" : "SHSMKYFFTS",
            "DRB1*01:01" : "RFLWQLKFEC",
            "DRB3*01:01" : "RFLELRKSEC",
            "DRB4*01:01" : "RFLEQAKCEC",
            "DRB5*01:01" : "RFLQQDKYEC",
            "DQA1*01:01" : "DHVASCGVNL",
            "DQB1*05:01" : "DFVYQFKGLC",
            "DPA1*01:03" : "DHVSTYAAFV",
            "DPB1*01:01" : "NYVYQGRQEC"
        }

        ard_first_ten_full = {
            "A*01:01:01:01" : "SHSMRYFFTS",
            "B*07:02:01:01" : "SHSMRYFYTS",
            "C*01:02:01:01" : "SHSMKYFFTS",
            "DRB1*01:01:01:01" : "RFLWQLKFEC",
            "DRB3*01:01:02:01" : "RFLELRKSEC",
            "DRB4*01:01:01:01" : "RFLEQAKCEC",
            "DRB5*01:01:01:01" : "RFLQQDKYEC",
            "DQA1*01:01:01:01" : "DHVASCGVNL",
            "DQB1*05:01:01:01" : "DFVYQFKGLC",
            "DPA1*01:03:01:01" : "DHVSTYAAFV",
            "DPB1*01:01:01:01" : "NYVYQGRQEC"
        }

        refseq_full = {
            "A" : "A*01:01:01:01",
            "B" : "B*07:02:01:01",
            "C" : "C*01:02:01:01",
            "DRB1" : "DRB1*01:01:01:01",
            "DRB3" : "DRB3*01:01:02:01",
            "DRB4" : "DRB4*01:01:01:01",
            "DRB5" : "DRB5*01:01:01:01",
            "DQA1" : "DQA1*01:01:01:01",
            "DQB1" : "DQB1*05:01:01:01",
            "DPA1" : "DPA1*01:03:01:01",
            "DPB1" : "DPB1*01:01:01:01",
        }

        refseq = {
            "A" : "A*01:01",
            "B" : "B*07:02",
            "C" : "C*01:02",
            "DRB1" : "DRB1*01:01",
            "DRB3" : "DRB3*01:01",
            "DRB4" : "DRB4*01:01",
            "DRB5" : "DRB5*01:01",
            "DQA1" : "DQA1*01:01",
            "DQB1" : "DQB1*05:01",
            "DPA1" : "DPA1*01:03",
            "DPB1" : "DPB1*01:01",
        }

        # use only ARD positions
        ard_start_pos = {
            "A" : 2,
            "B" : 2,
            "C" : 2,
            "DRB1" : 6,
            "DRB3" : 6,
            "DRB4" : 6,
            "DRB5" : 6,
            "DQA1" : 6,
            "DQB1" : 6,
            "DPA1" : 4,
            "DPB1" : 6,
        }
        ard_end_pos = {
            "A" : 182,
            "B" : 182,
            "C" : 182,
            "DRB1" : 94,
            "DRB3" : 94,
            "DRB4" : 94,
            "DRB5" : 94,
            "DQA1" : 87,
            "DQB1" : 94,
            "DPA1" : 84,
            "DPB1" : 92,
        }

        # use cDNAcoordinates in hla.xml to compute
        # HLA-A - Exons 2, 3, 4 are coordinates 74 to 895
        # (895-74+1) / 3 = 274
        xrd_start_pos = {
            "A" : 2,
            "B" : 2,
            "C" : 2,
            "DRB1" : 6,
            "DRB3" : 6,
            "DRB4" : 6,
            "DRB5" : 6,
            "DQA1" : 6,
            "DQB1" : 6,
            "DPA1" : 4,
            "DPB1" : 6,
        }

        xrd_end_pos = {
            "A" : 274,
            "B" : 274,
            "C" : 274,
            "DRB1" : 188,
            "DRB3" : 188,
            "DRB4" : 188,
            "DRB5" : 188,
            "DQA1" : 181,
            "DQB1" : 188,
            "DPA1" : 178,
            "DPB1" : 186,
        }



        self.refseq_full = refseq_full
        self.refseq = refseq
        self.ard_start_pos = ard_start_pos
        self.ard_end_pos = ard_end_pos
        self.xrd_start_pos = xrd_start_pos
        self.xrd_end_pos = xrd_end_pos
        self.hlaProteinOffset = hlaProteinOffset
        self.first_ten = first_ten_full
        self.main()

        if self.ungap == True:
            self.remove_gap()
            self.ard_end_pos = ard_end_pos

    # generate a regular expression to find first ten residues of each locus
    # ignoring '-' characters.
    def regex_gen(self):
        regexes = {}
        for each in self.first_ten.keys():
            regex = ""
            for i in range(0,10):
                regex += str(self.first_ten[each][i])+"[^"+self.first_ten[each]+"]*?"
                reg_loc = each.split('*')[0]
            regexes[reg_loc] = regex
        return regexes

    # Align the coordinate system as appropriate to the beginning of the mature
    # protein sequence.
    def coordinate(self, regex, sequence):
        o = re.search(regex, str(sequence))
        offset = o.start()
        return offset

    def getMatureProteinOffset(self, locus):
            return self.hlaProteinOffset.get(locus, "Invalid HLA Locus")

    def adjust_end(self, multipleseq,loc,ard_start_pos,ard_end_pos,prev=0):
        loc_full_allele = self.refseq_full[loc]
        full_protein = ''
        while full_protein == '':
            full_protein = multipleseq[id==loc_full_allele].seq
            if full_protein == '':
                loc_full_allele = ':'.join(loc_full_allele.split(':')[:-1])
        self.refseq_full[loc] = loc_full_allele
        mature_protein = full_protein[self.getMatureProteinOffset(loc):]
        start = ard_start_pos
        end = ard_end_pos
        count = mature_protein[start:end].count('-')

        check = count - prev

        if (check == 0):
            new_end = end + (count-prev)
            return new_end
        else:
            new_end = end + (count - prev)
            prev = count
            return self.adjust_end(multipleseq,loc,start,new_end,prev)

    def remove_ins(self, loc_full_alseq):
        # need to remove the inserts from the reference sequence to print into
        # the IMGT/HLA .txt file
        gapframe = pd.DataFrame.from_dict(loc_full_alseq, orient="index")
        droplist = []
        for i, row in gapframe.iterrows():
            if i == self.refseq[loc]:
                for name, data in gapframe.iteritems():
                    if data[i] == '-':
                        droplist.append(name)
                    else:
                        continue
            else:
                continue
        ungapframe = gapframe.drop(droplist, axis=1)
        for j, jrow in ungapframe.iterrows():
            jrow = jrow.to_string(header=False, index=False)
            jrow = jrow.replace('\n', '')
            jrow = jrow.replace(' ', '')
            loc_full_alseq[j] = jrow

        return loc_full_alseq

    def reference(self, multipleseq, loc):
        loc_full_allele = self.refseq_full[loc]
        regexes = self.regex_gen()
        mslist = []
        while mslist == []:
            mslist = [x for x in multipleseq if x.id == loc_full_allele]
            if mslist == []:
                loc_full_allele = ':'.join(loc_full_allele.split(':')[:-1])
        self.refseq_full[loc] = loc_full_allele
        loc_full_protein = mslist[0]
        #loc_full_protein = multipleseq[id==loc_full_allele]
        full_protein = loc_full_protein.seq

        # align coordinate system based on reference sequence.
        # uses two field allele if full allele not present.
        # does not overwrite full allele if present.
        # necessary for earlier database versions.
        offset = self.coordinate(regexes[loc], full_protein)
        self.hlaProteinOffset[loc] = offset

        return

    def adjust(self, multipleseq, loc):
        # then adjust end position based on dashes.
        new_end = self.adjust_end(multipleseq, loc, self.ard_start_pos[loc],
                                            self.ard_end_pos[loc],prev=0)
        self.ard_end_pos[loc] = new_end
        new_xrd_end = self.adjust_end(multipleseq, loc, self.xrd_start_pos[loc],
                                            self.xrd_end_pos[loc],prev=0)
        self.xrd_end_pos[loc] = new_xrd_end
        return


    def generate_IMGT(self, HLA_full_alseq):
        outfile = open("./IMGT_HLA_Full_Protein_" + str(self.dbversion) + ".txt", "w+")
        outfile.write("Allele\tFull_Protein\n")
        for allele_loctype in HLA_full_alseq:
            outfile.write("HLA-" + allele_loctype + "\t" +
                        str(HLA_full_alseq[allele_loctype]) + "\n")
        outfile.close()
        return

    def remove_gap(self):
        loci =  ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        refgaps2 = {}
        refgaps4 = {}
        refgaps3 = {}
        for locus in loci:
            seq2 = str(self.HLA_seq[self.refseq[locus]].seq)
            seq4 = str(self.HLA_full_allele[self.refseq_full[locus]].seq)
            seq3 = str(self.HLA_all_records[self.refseq_full[locus]].seq)
            refgaps2[locus] = [i for i, ltr in enumerate(seq2) if ltr == '-']
            refgaps4[locus] = [j for j, lttr in enumerate(seq4) if lttr == '-']
            refgaps3[locus] = [k for k, ltttr in enumerate(seq3) if ltttr == '-']
        for x2 in self.HLA_seq.keys():
            locus = x2.split('*')[0]
            gapseq2 = str(self.HLA_seq[x2].seq)
            seq2x = ''.join([gapseq2[x] for x in range(len(gapseq2)) if not((gapseq2[x] == '-') and (x in refgaps2[locus]))])
            self.HLA_seq[x2].seq = Seq(seq2x)
        for x4 in self.HLA_full_allele.keys():
            locus = x4.split('*')[0]
            gapseq4 = str(self.HLA_full_allele[x4].seq)
            seq4x = ''.join([gapseq4[y] for y in range(len(gapseq4)) if not((gapseq4[y] == '-') and (y in refgaps4[locus]))])
            self.HLA_full_allele[x4].seq = Seq(seq4x)
        # added this loop for validation purposes
        for x3 in self.HLA_all_records.keys():
            locus = x3.split('*')[0]
            gapseq3 = str(self.HLA_all_records[x3].seq)
            seq3x = ''.join([gapseq3[z] for z in range(len(gapseq3)) if not((gapseq3[z] == '-') and (z in refgaps3[locus]))])
            self.HLA_all_records[x3].seq = Seq(seq3x)
        return

    def main(self):
        loci =  ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        if self.imputed == True:
            loci =  ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        self.HLA_full_allele = {} # Full four-field allele names
        #self.HLA_full_alseq = {} # Only used for generating IMGTHLA files
        self.HLA_seq = {} # Two-field
        self.HLA_all_records = {} # for validation
        regex = '\w*\*\d*\:\d*'
        suffixes = ["L", "S", "C", "A", "Q", "N"]
        fdir = self.fdir
        ftype = self.ftype
        for locus in loci:
            loc_full_alseq = {}
            if ftype == "msf":
                seq_filename = "{}/{}/{}_prot_{}.{}".format(str(pathloc), fdir, locus, str(self.dbversion), ftype)
                seq_filename = Path(seq_filename)
            elif ftype == "fasta":
                seq_filename = "{}/{}/{}_{}.{}".format(str(pathloc), fdir, locus, str(self.dbversion), ftype)
                seq_filename = Path(seq_filename)
            if (path.exists(seq_filename) == False) & (ftype == 'msf'):
                print("Downloading requested MSF files for locus " + locus + "...")
                url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/" + str(
                    self.dbversion) + "/msf/" + locus + "_prot.msf"
                r = requests.get(url)
                with open(seq_filename, 'wb') as f:
                    f.write(r.content)
                multipleseq = AlignIO.read(seq_filename, format="msf")
            elif (path.exists(seq_filename) == False) & (ftype == 'fasta'):
                print('ERROR - NO IMPUTED SEQUENCE FILES PRESENT - RUN IMPUTATION FIRST')
                break
            else:
                if ftype == 'msf':
                    multipleseq = AlignIO.read(seq_filename, format="msf")
                elif ftype == 'fasta':
                    print(locus)
                    multipleseq = SeqIO.parse(seq_filename, format='fasta')
                    multipleseq = list(multipleseq)
                    # temporarily removing alleles with indels
                    multipleseq = [x for x in multipleseq if x.id[-1] != 'Q']
                    print("MAX:\t"+str(max(len(sequ.seq) for sequ in multipleseq)))
                    print("MIN:\t"+str(min(len(sequ.seq) for sequ in multipleseq)))
                    #multipleseq = AlignIO.read(seq_filename, format="fasta")

            #print("OLD {} XRD END:\t{}".format(locus,str(self.xrd_end_pos[locus])))

            if locus == 'DRB345':
                for loc in ['DRB3', 'DRB4', 'DRB5']:
                    self.reference(multipleseq, loc)
                    self.adjust(multipleseq, loc)
            else:
                self.reference(multipleseq, locus)
                self.adjust(multipleseq, locus)

            #print("NEW {} XRD END:\t{}".format(locus,str(self.xrd_end_pos[locus])))

            for record in multipleseq:
                loc_full_allele = record.id
                (loc,full_allele) = loc_full_allele.split("*")
                # append the suffix - needed for null alleles
                # index starts at 1 since some loci characters are also suffixes
                separator = loc_full_allele.find("*")
                # TEMP - not allowing for two-field nulls
                if any(x in loc_full_allele[separator:] for x in suffixes):
                    #loc_two_field_allele = re.match(regex, loc_full_allele).group() + \
                    #                    loc_full_allele[-1]
                    next
                else:
                    try:
                        loc_two_field_allele = re.match(regex, loc_full_allele).group()
                    except:
                        print(loc_full_allele)
                full_protein = record.seq
                nogap = full_protein

                loc_full_alseq[loc_full_allele] = SeqRecord(nogap)

                # skip missing sequences in IMGT/HLA file
                if (len(full_protein) <10): # #NAME?
                    print ("Missing Sequence:" + loc_full_allele)
                    continue

                mature_protein = full_protein[self.getMatureProteinOffset(loc):]
                mrecord = SeqRecord(mature_protein)


                # full allele name
                self.HLA_full_allele[loc_full_allele] = mrecord

                # don't overwrite two-field alleles with new sequences - more likely
                # to be incomplete
                if (loc_two_field_allele not in self.HLA_seq.keys()):
                    self.HLA_seq[loc_two_field_allele] = mrecord
                # print (HLA_seqrecord_dict[allele])

                # include all records for validation purposes
                self.HLA_all_records[loc_full_allele] = mrecord
                if loc_full_allele.count(':') == 3:
                    loc_three_field = ':'.join(loc_full_allele.split(':')[:-1])
                if loc_three_field not in self.HLA_all_records.keys():
                    self.HLA_all_records[loc_three_field] = mrecord

                # TODO - add feature annotation to SeqIO object
                # https://biopython.org/wiki/SeqRecord
                # e.g. - which allele contain Bw4/Bw6 epitopes, bind LILRB1, etc
                # https://www.biostars.org/p/57549/

                # print (HLA_seqrecord_dict[allele].seq)

            #!GB!# Commented out these two lines as well as the generate_IMGT()
            # function call, since the process is intensive and does not need to be
            #!GB!# run every time the code is used. Might should consider a switch to
            # splitting a HLA_full_alseq DataFrame by locus and doing all of the work in
            #!GB!# a single function definition.
            #loc_full_alseq = remove_ins(loc_full_alseq)
            #HLA_full_alseq.update(loc_full_alseq)

        #generate_IMGT(HLA_full_alseq)

        return

    def getAAsubstring(self, allele,start_position,end_position):
        try:
            AAsubstring = self.HLA_seq[allele].seq[start_position-1:end_position]
        except KeyError:
            AAsubstring = self.HLA_full_allele[allele].seq[start_position-1:end_position]
        return AAsubstring

    # get the AA at any position of any HLA allele
    def getAAposition(self, allele,position):
        position -= 1
        try:
            AAposition = self.HLA_seq[allele].seq[position]
        except KeyError:
            AAposition = self.HLA_full_allele[allele].seq[position]
        return AAposition


    # print amino acid mismatches between two alleles
    def AllelePair_ARD_AA_Matched(self, allele1, allele2):
        (loc,_allele) = allele1.split('*')
        start = self.ard_start_pos[loc]
        end = self.ard_end_pos[loc]
        MM_any = 0

        # position-by-position is too slow
        # for position in range(start,end):
        #     MM = self.isPositionMismatched(allele1,allele2,position)
        #     if (MM == 1):
        #         MM_any = 1

        # compare entire two-field ARD string
        # if (self.HLA_seq[allele1].seq != self.HLA_seq[allele2].seq):
        #     MM_any = 1

        # TODO - Use sequences with missing data filled in
        allele1_ARDseq = self.getAAsubstring(allele1,start,end)
        allele2_ARDseq = self.getAAsubstring(allele2,start,end)
        # print (self.HLA_seq[allele1].seq)
        # print (self.HLA_seq[allele2].seq)
        # print (allele1_ARDseq)
        # print (allele2_ARDseq)
        # print (str(len(allele1_ARDseq)))
        if (allele1_ARDseq != allele2_ARDseq):
            MM_any = 1

        return (MM_any)

    # print amino acid mismatches between two alleles (XRD)
    def AllelePair_XRD_AA_Matched(self, allele1, allele2):
        (loc,_allele) = allele1.split('*')
        start = self.xrd_start_pos[loc]
        end = self.xrd_end_pos[loc]
        MM_any = 0

        # position-by-position is too slow
        # for position in range(start,end):
        #     MM = self.isPositionMismatched(allele1,allele2,position)
        #     if (MM == 1):
        #         MM_any = 1

        # compare entire two-field ARD string
        # if (self.HLA_seq[allele1].seq != self.HLA_seq[allele2].seq):
        #     MM_any = 1

        # TODO - Use sequences with missing data filled in
        allele1_XRDseq = self.getAAsubstring(allele1,start,end)
        allele2_XRDseq = self.getAAsubstring(allele2,start,end)
        # print (self.HLA_seq[allele1].seq)
        # print (self.HLA_seq[allele2].seq)
        # print (allele1_ARDseq)
        # print (allele2_ARDseq)
        # print (str(len(allele1_ARDseq)))
        if (allele1_XRDseq != allele2_XRDseq):
            MM_any = 1

        return (MM_any)

    # get the SFVT epitope name (with underscore)
    def getEpitope(self, allele,position_list):
        sfvt_list = []
        for position in position_list:
            AA = self.getAAposition(allele,position)
            sfvt_aa = str(position) + AA
            sfvt_list.append(sfvt_aa)

        sfvt = "_".join(sfvt_list)
        return (sfvt)

    # given two alleles at same locus and a position - Yes/No if mismatched
    def isPositionMismatched(self, allele1,allele2,position):
        AA_allele1 = self.getAAposition(allele1,position)
        AA_allele2 = self.getAAposition(allele2,position)
        if (AA_allele1 != AA_allele2):
            return True
        else:
            return False

    # count number of mismatches at position between donor and recip
    # getTruth function from runMatchMC
    def count_AA_Mismatches(self, aa1_donor,aa2_donor,aa1_recip,aa2_recip):
        mm_count = 0
        if ((aa1_donor != aa1_recip) & (aa1_donor != aa2_recip)):
            mm_count+=1
        if ((aa2_donor != aa1_recip) & (aa2_donor != aa2_recip)):
            mm_count+=1
        return mm_count

    # count number of mismatches between alleles at a given position, adjusting
    # for donor homozygosity
    def count_AA_Mismatches_Allele(self, allele1_donor,allele2_donor,allele1_recip,
                                allele2_recip,position):
        donor_homoz = 0
        if (allele1_donor == allele2_donor):
            donor_homoz = 1
        aa1_donor = self.getAAposition(allele1_donor,position)
        aa2_donor = self.getAAposition(allele2_donor,position)
        aa1_recip = self.getAAposition(allele1_recip,position)
        aa2_recip = self.getAAposition(allele2_recip,position)

        mm_count = self.count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

        if ((mm_count == 2) & (donor_homoz == 1)):
            mm_count = 1

        return mm_count

    def count_AA_Mismatches_SFVT(self, allele1_donor,allele2_donor,allele1_recip,
                                allele2_recip,position_list):
        donor_homoz = 0
        mm_total = 0
        if (allele1_donor == allele2_donor):
            donor_homoz = 1

        # increment MM count for all positions in list
        for position in position_list:
            aa1_donor = self.getAAposition(allele1_donor,position)
            aa2_donor = self.getAAposition(allele2_donor,position)
            aa1_recip = self.getAAposition(allele1_recip,position)
            aa2_recip = self.getAAposition(allele2_recip,position)

            mm_count = self.count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

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
    def AA_MM_Allele(self, allele1_donor,allele2_donor,allele1_recip,allele2_recip,
                    position):
        aa1_donor = self.getAAposition(allele1_donor,position)
        aa2_donor = self.getAAposition(allele2_donor,position)
        aa1_recip = self.getAAposition(allele1_recip,position)
        aa2_recip = self.getAAposition(allele2_recip,position)

        is_mm = self.AA_MM(aa1_donor,aa2_donor,aa1_recip,aa2_recip)

        return is_mm

    # weighted choice from https://scaron.info/blog/python-weighted-choice.html
    def weighted_choice(self, seq, weights):
        assert len(weights) == len(seq)
        assert abs(1. - sum(weights)) < 1e-6

        x = random.random()
        for i, elmt in enumerate(seq):
            if x <= weights[i]:
                return elmt
            x -= weights[i]
