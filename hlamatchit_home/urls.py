from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name = 'hlamatchit_home'),
    path('about/', views.about, name = 'hlamatchit_home/about'),
    path('allele2aa/', views.allele2aa, name = 'hlamatchit_home/allele2aa'),
]
