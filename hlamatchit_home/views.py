from urllib import response
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched


# Create your views here.

def allele2aa(request):
    return render(request, 'hlamatchit_home/allele2aa.html')

def allele2aa_out(request):
    allele = request.GET['userinput1']
    position = request.GET['userinput2']

    print (allele + " " + position)

    AA = getAAposition(allele,position)

    return render(request, 'hlamatchit_home/allele2aa_out.html', 
        {'Allele': allele, 'Position': position, 'AA': AA})

def ismatch(request):
    return render(request, 'hlamatchit_home/ismatch.html')

def ismatch_out(request):
    allele1 = request.GET['userinput1']
    allele2 = request.GET['userinput2']
    position = request.GET['userinput3']

    print (allele1 + " " + allele2 + " " + position)

    Match = isPositionMismatched(allele1, allele2,position)

    return render(request, 'hlamatchit_home/ismatch_out.html', 
        {'Allele 1': allele1,'Allele 2': allele2, 'Position': position,'Position is Matched': Match })


def home(request):
    return render(request, 'hlamatchit_home/home.html')

def about(request):
    return render(request, 'hlamatchit_home/about.html')


