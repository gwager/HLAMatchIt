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
    race1 = request.GET['userinput17']
    antigen1 = request.GET['userinput18']
    race2 = request.GET['userinput19']
    antigen2 = request.GET['userinput20']

    allele1 = antigen2HFallele(race1,antigen1)
    print(allele1)

    allele2 = antigen2HFallele(race2,antigen2)
    print(allele2)

    return render(request, 'antigen2aa_out.html', 
        {'Race1': race1,'Antigen1': antigen1,'Allele1': allele1,'Race2': race2,'Antigen2': antigen2,'Allele2': allele2})


