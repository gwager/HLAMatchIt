import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_fibers import antigen2HFallele,getAAgenostringmatch, afibershazard, cfibershazard, bfibershazard, drfibershazard, dqfibershazard


# Create your views here.

def antigen2aa(request):
    return render(request, 'antigen2aa.html')

def antigen2aa_out(request):
    drace = request.GET['userinput17']
    danta1 = request.GET['userinput18']
    danta2 = request.GET['userinput23']
    dantc1 = request.GET['userinput19']
    dantc2 = request.GET['userinput24']
    dantb1 = request.GET['userinput20']
    dantb2 = request.GET['userinput25']
    dantdr1 = request.GET['userinput21']
    dantdr2 = request.GET['userinput26']
    dantdq1 = request.GET['userinput22']
    dantdq2 = request.GET['userinput27']

    rrace = request.GET['userinput28']
    ranta1 = request.GET['userinput29']
    ranta2 = request.GET['userinput34']
    rantc1 = request.GET['userinput30']
    rantc2 = request.GET['userinput35']
    rantb1 = request.GET['userinput31']
    rantb2 = request.GET['userinput36']
    rantdr1 = request.GET['userinput32']
    rantdr2 = request.GET['userinput37']
    rantdq1 = request.GET['userinput33']
    rantdq2 = request.GET['userinput38']


    aloc = "A"
    cloc = "C"
    bloc = "B"
    drloc = "DRB1"
    dqloc = "DQB1"

    #get Donor most prob alleles
    dalla1,dfa1 = antigen2HFallele(drace,danta1)
    dalla2,dfa2 = antigen2HFallele(drace,danta2)

    dallc1,dfc1 = antigen2HFallele(drace,dantc1)
    dallc2,dfc2 = antigen2HFallele(drace,dantc2)

    dallb1,dfb1 = antigen2HFallele(drace,dantb1)
    dallb2,dfb2 = antigen2HFallele(drace,dantb2)

    dalldr1,dfdr1 = antigen2HFallele(drace,dantdr1)
    dalldr2,dfdr2 = antigen2HFallele(drace,dantdr2)

    dalldq1,dfdq1 = antigen2HFallele(drace,dantdq1)
    dalldq2,dfdq2 = antigen2HFallele(drace,dantdq2)


    #get Recipient most prob alleles
    ralla1,rfa1 = antigen2HFallele(rrace,ranta1)
    ralla2,rfa2 = antigen2HFallele(rrace,ranta2)

    rallc1,rfc1 = antigen2HFallele(rrace,rantc1)
    rallc2,rfc2 = antigen2HFallele(rrace,rantc2)

    rallb1,rfb1 = antigen2HFallele(rrace,rantb1)
    rallb2,rfb2 = antigen2HFallele(rrace,rantb2)

    ralldr1,rfdr1 = antigen2HFallele(rrace,rantdr1)
    ralldr2,rfdr2 = antigen2HFallele(rrace,rantdr2)

    ralldq1,rfdq1 = antigen2HFallele(rrace,rantdq1)
    ralldq2,rfdq2 = antigen2HFallele(rrace,rantdq2)

    dAa1,dAa2,rAa1,rAa2,Acount,Apos,Apos1,Aend_pos = getAAgenostringmatch(dalla1,dalla2,ralla1,ralla2, aloc)
    dCa1,dCa2,rCa1,rCa2,Ccount,Cpos,Cpos1,Cend_pos = getAAgenostringmatch(dallc1,dallc2,rallc1,rallc2, cloc)
    dBa1,dBa2,rBa1,rBa2,Bcount,Bpos,Bpos1,Bend_pos = getAAgenostringmatch(dallb1,dallb2,rallb1,rallb2, bloc)
    dDRa1,dDRa2,rDRa1,rDRa2,DRcount,DRpos,DRpos1,DRend_pos = getAAgenostringmatch(dalldr1,dalldr2,ralldr1,ralldr2, drloc)
    dDQa1,dDQa2,rDQa1,rDQa2,DQcount,DQpos,DQpos1,DQend_pos = getAAgenostringmatch(dalldq1,dalldq2,ralldq1,ralldq2, dqloc)


    arange = Aend_pos
    crange = Cend_pos
    brange = Bend_pos
    drrange = DRend_pos
    dqrange = DQend_pos

    ahazard, Apos2 = afibershazard(Apos)
    chazard, Cpos2 = cfibershazard(Cpos)
    bhazard, Bpos2 = bfibershazard(Bpos)
    drhazard, DRpos2, drprob = drfibershazard(DRpos)
    dqhazard, DQpos2, dqprob = dqfibershazard(DQpos)
    #fprob = fibersprobability(locus, pos)

    
    return render(request, 'antigen2aa_out.html', 
        {'drace': drace,'danta1': danta1,'dalla1': dalla1, 'dfa1': dfa1,'danta2': danta2,'dalla2': dalla2,'dfa2': dfa2, 'dantc1': dantc1, 'dallc1': dallc1 ,'dfc1': dfc1, 'dantc2':dantc2,'dallc2': dallc2, 'dfc2': dfc2,'dantb1': dantb1,'dallb1': dallb1, 'dfb1': dfb1, 'dantb2': dantb2,'dallb2': dallb2, 'dfb2': dfb2,'dantdr1': dantdr1, 'dalldr1': dalldr1, 'dfdr1': dfdr1, 'dantdr2': dantdr2, 'dalldr2': dalldr2, 'dfdr2' : dfdr2, 'dantdq1': dantdq1, 'dalldq1':dalldq1,'dfdq1':dfdq1,'dantdq2':dantdq2,'dalldq2':dalldq2,'dfdq2':dfdq2, 'rrace': rrace,'ranta1': ranta1,'ralla1': ralla1, 'rfa1': rfa1,'ranta2': ranta2,'ralla2': ralla2,'rfa2': rfa2, 'rantc1': rantc1, 'rallc1': rallc1 ,'rfc1': rfc1, 'rantc2':rantc2,'rallc2': rallc2, 'rfc2': rfc2,'rantb1': rantb1,'rallb1': rallb1, 'rfb1': rfb1, 'rantb2': rantb2,'rallb2': rallb2, 'rfb2': rfb2,'rantdr1': rantdr1, 'ralldr1': ralldr1, 'rfdr1': rfdr1, 'rantdr2': rantdr2, 'ralldr2': ralldr2, 'rfdr2' : rfdr2,'rantdq1': rantdq1, 'ralldq1':ralldq1,'rfdq1':rfdq1,'rantdq2':rantdq2,'ralldq2':ralldq2,'rfdq2':rfdq2, 'arange': arange, 'Acount':Acount, 'Apos' : Apos, 'Apos1': Apos1, 'Apos2' : Apos2, 'ahazard' : ahazard, 'crange': crange, 'Ccount':Ccount, 'Cpos' : Cpos, 'Cpos1': Cpos1, 'Cpos2' : Cpos2, 'chazard' : chazard,'brange': brange, 'Bcount':Bcount, 'Bpos' : Bpos, 'Bpos1': Bpos1, 'Bpos2' : Bpos2, 'bhazard' : bhazard,'drrange': drrange, 'DRcount':DRcount, 'DRpos' : DRpos, 'DRpos1': DRpos1, 'DRpos2' : DRpos2, 'drhazard' : drhazard, 'drprob' : drprob, 'dqrange': dqrange, 'DQcount':DQcount, 'DQpos' : DQpos, 'DQpos1': DQpos1, 'DQpos2' : DQpos2, 'dqhazard' : dqhazard, 'dqprob' : dqprob})
