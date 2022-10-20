import sys
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched


# Create your views here.

def allele2aa(request):
    return render(request, 'allele2aa.html')

def allele2aa_out(request):
    allele = request.GET['userinput1']
    position = request.GET['userinput2']
    
    print (allele + " " + position)

    AA = getAAposition(allele,position)

    return render(request, 'allele2aa_out.html', 
        {'Allele': allele, 'Position': position, 'AA': AA})


