from django.urls import path
from django.conf import settings
from django.conf.urls import include
from django.contrib import admin
from urllib import response
from django.shortcuts import render
from hlamatchit_home import views
from . import views
from . import views_aa
from . import views_string
from . import views_match
from . import views_mmcount
from . import views_stringcount
from . import views_antigen2aa
from . import views_DRQantigen2aa



urlpatterns = [
    path('', views.home, name = 'hlamatchit_home'),
    path('about/', views.about, name = 'about.html'),
    path('allele2aa.html/', views_aa.allele2aa, name = 'allele2aa.html'),
    path('allele2aa_out.html/', views_aa.allele2aa_out, name = 'allele2aa_out.html'),
    path('allele2string.html/', views_string.allele2string, name = 'allele2string.html'),
    path('allele2string_out.html/', views_string.allele2string_out, name = 'allele2string_out.html'),
    path('ismatch.html/', views_match.ismatch, name = 'ismatch.html'),
    path('ismatch_out.html/', views_match.ismatch_out, name = 'ismatch_out.html'),
    path('mmcount.html/', views_mmcount.mmcount, name = 'mmcount.html'),
    path('mmcount_out.html/', views_mmcount.mmcount_out, name = 'mmcount_out.html'),
    path('stringcount.html/', views_stringcount.stringcount, name = 'stringcount.html'),
    path('stringcount_out.html/', views_stringcount.stringcount_out, name = 'stringcount_out.html'),
    path('antigen2aa.html/', views_antigen2aa.antigen2aa, name = 'antigen2aa.html'),
    path('antigen2aa_out.html/', views_antigen2aa.antigen2aa_out, name = 'antigen2aa_out.html'),
    path('DRQantigen2aa.html/', views_DRQantigen2aa.DRQantigen2aa, name = 'DRQantigen2aa.html'),
    path('DRQantigen2aa_out.html/', views_DRQantigen2aa.DRQantigen2aa_out, name = 'DRQantigen2aa_out.html'),
]

