import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched
from .aa_fibers import getAAstringmatch


# Create your views here.

def stringcount(request):
    return render(request, 'stringcount.html')

def stringcount_out(request):
    aa1_donor = request.GET['userinput13']
    aa2_donor = request.GET['userinput14']
    start_pos = request.GET['userinput15']
    end_pos = request.GET['userinput16']

    range = end_pos

    a1,a2,count,pos = getAAstringmatch(aa1_donor,aa2_donor,start_pos,end_pos)
    print(count)

    return render(request, 'stringcount_out.html', 
        {'Type1': aa1_donor,'Type2': aa2_donor,'Range': range ,'str1': a1, 'str2':a2,'total': count, 'pos': pos })


