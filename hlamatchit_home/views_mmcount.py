import sys
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched
from .aa_matching import count_AA_Mismatches


# Create your views here.

def mmcount(request):
    return render(request, 'mmcount.html')

def mmcount_out(request):
    aa1_donor = request.GET['userinput6']
    aa2_donor = request.GET['userinput7']
    aa1_recip = request.GET['userinput8']
    aa2_recip = request.GET['userinput9']

    

    count = count_AA_Mismatches(aa1_donor,aa2_donor,aa1_recip,aa2_recip)
    print(count)

    return render(request, 'mmcount_out.html', 
        {'TypeA1': aa1_donor,'TypeA2': aa2_donor, 'TypeB1': aa1_recip,'TypeB2': aa2_recip, 'Mismatches': count })


