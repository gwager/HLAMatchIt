from django.shortcuts import render
from django.http import HttpResponse
# Create your views here.

def home(request):
    return render(request, 'hlamatchit_home/home.html')

def about(request):
    return render(request, 'hlamatchit_home/about.html')

def allele2aa(request):
    return render(request, 'hlamatchit_home/allele2aa.html')

