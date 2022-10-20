from django.urls import path
from django.conf import settings
from django.conf.urls import include
from django.contrib import admin
from urllib import response
from django.shortcuts import render
from . import views
from . import views_aa
from . import views_match


urlpatterns = [
    path('', views.home, name = 'hlamatchit_home'),
    path('about/', views.about, name = 'about.html'),
    path('allele2aa.html/', views_aa.allele2aa, name = 'allele2aa.html'),
    path('allele2aa_out.html/', views_aa.allele2aa_out, name = 'allele2aa_out.html'),
    path('ismatch.html/', views_match.ismatch, name = 'ismatch.html'),
    path('ismatch_out.html/', views_match.ismatch_out, name = 'ismatch_out.html'),
]

