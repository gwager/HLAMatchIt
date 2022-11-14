import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_fibers import getAAsubstring


# Create your views here.

def allele2string(request):
    return render(request, 'allele2string.html')

def allele2string_out(request):
    allele = request.GET['userinput10']
    start = request.GET['userinput11']
    end = request.GET['userinput12']

    

    string = str(getAAsubstring(allele,start,end))

    print(string)

    return render(request, 'allele2string_out.html', 
        {'Allele': allele, 'Start': start, 'End': end ,'aastring': string})


