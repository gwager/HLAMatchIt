import sys
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched


# Create your views here.

def home(request):
    return render(request, 'home.html')

def about(request):
    return render(request, 'about.html')


