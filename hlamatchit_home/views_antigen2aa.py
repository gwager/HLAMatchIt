import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched
from .aa_fibers import antigen2HFallele, antigen2allele, highfreq, freqfileselect,getAAstringmatch


# Create your views here.

def antigen2aa(request):
    return render(request, 'antigen2aa.html')

def antigen2aa_out(request):
    locus = request.GET['userinput17']
    race1 = request.GET['userinput18']
    antigen1 = request.GET['userinput19']
    race2 = request.GET['userinput20']
    antigen2 = request.GET['userinput21']


    allele1,freq1 = antigen2HFallele(race1,antigen1)
    #print(allele1,freq1)

    print(antigen2)
    allele2,freq2 = antigen2HFallele(race2,antigen2)
    #print(allele2)
    a1,a2,count,pos,pos1,pos2, end_pos = getAAstringmatch(allele1,allele2,locus)

    range = end_pos

    
    return render(request, 'antigen2aa_out.html', 
        {'Race1': race1,'Antigen1': antigen1,'Allele1': allele1, 'Freq1': freq1,'Race2': race2,'Antigen2': antigen2,'Allele2': allele2, 'Freq2': freq2, 'Range': range ,'str1': a1, 'str2':a2,'total': count, 'pos': pos,'pos1': pos1,'pos2': pos2 })


