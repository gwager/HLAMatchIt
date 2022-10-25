import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched
from .aa_fibers import antigen2HFallele, antigen2allele, highfreq, freqfileselect,getAAstringmatch, AA_MM


# Create your views here.

def antigen2aa(request):
    return render(request, 'antigen2aa.html')

def antigen2aa_out(request):
    race1 = request.GET['userinput17']
    antigen1 = request.GET['userinput18']
    race2 = request.GET['userinput19']
    antigen2 = request.GET['userinput20']

    allele1,freq1 = antigen2HFallele(race1,antigen1)
    #print(allele1,freq1)

    print(antigen2)
    allele2,freq2 = antigen2HFallele(race2,antigen2)
    #print(allele2)

    loc = ['A','C','B','DRB1','DQA1','DQB1','DPA1','DPB1']
    for loc in allele1:
        if (loc == 'A'):
        Aloc = GENO_A[0].split('*')[0]
        start = aa_mm.ard_start_pos[Aloc]
        end = aa_mm.ard_end_pos[Aloc]

    return render(request, 'antigen2aa_out.html', 
        {'Race1': race1,'Antigen1': antigen1,'Allele1': allele1, 'Freq1': freq1,'Race2': race2,'Antigen2': antigen2,'Allele2': allele2, 'Freq2': freq2})


