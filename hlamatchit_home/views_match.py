import sys
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_fibers import isPositionMismatched


# Create your views here.

def ismatch(request):
    return render(request, 'ismatch.html')

def ismatch_out(request):
    allele1 = request.GET['userinput3']
    allele2 = request.GET['userinput4']
    position = request.GET['userinput5']

    

    Match = isPositionMismatched(allele1,allele2,position)
    print(Match)

    if(Match == 0):
        Match = "TRUE"
    else:
        Match = "False"
    print(Match)
    print (allele1 + " " + allele2 + " " + position)

    return render(request, 'ismatch_out.html', 
        {'Allele1': allele1,'Allele2': allele2, 'Position': position,'Match': Match })


