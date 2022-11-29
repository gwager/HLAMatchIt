import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_fibers import antigen2HFalleleFIBERS, getAAgenostringmatchFIBERS,drfibersprob, dqfibersprob, getAAgenostringmatch


# Create your views here.

def DRQantigen2aa(request):
    return render(request, 'DRQantigen2aa.html')

def DRQantigen2aa_out(request):
    drace = request.GET['userinput39']

    dantdr1 = request.GET['userinput40']
    dantdr2 = request.GET['userinput42']
    dantdq1 = request.GET['userinput41']
    dantdq2 = request.GET['userinput43']

    rrace = request.GET['userinput44']
    
    rantdr1 = request.GET['userinput45']
    rantdr2 = request.GET['userinput47']
    rantdq1 = request.GET['userinput46']
    rantdq2 = request.GET['userinput48']


    drloc = "DRB1"
    dqloc = "DQB1"

    #get Donor most prob alleles

    dalldr1,dfdr1,dnalldr1, ddr1sumant, ddr1alleles= antigen2HFalleleFIBERS(drace,dantdr1)
    dalldr2,dfdr2,dnalldr2, ddr2sumant, ddr2alleles= antigen2HFalleleFIBERS(drace,dantdr2)

    dalldq1,dfdq1,dnalldq1, ddq1sumant, ddq1alleles= antigen2HFalleleFIBERS(drace,dantdq1)
    dalldq2,dfdq2,dnalldq2, ddq2sumant, ddq2alleles= antigen2HFalleleFIBERS(drace,dantdq2)


    #get Recipient most prob alleles

    ralldr1,rfdr1,rnalldr1, rdr1sumant, rdr1alleles= antigen2HFalleleFIBERS(rrace,rantdr1)
    ralldr2,rfdr2,rnalldr2, rdr2sumant, rdr2alleles = antigen2HFalleleFIBERS(rrace,rantdr2)

    ralldq1,rfdq1,rnalldq1, rdq1sumant,rdq1alleles = antigen2HFalleleFIBERS(rrace,rantdq1)
    ralldq2,rfdq2,rnalldq2, rdq2sumant, rdq2alleles= antigen2HFalleleFIBERS(rrace,rantdq2)

    print("Got alleles")

    drfprob = getAAgenostringmatchFIBERS(ddr1alleles,ddr2alleles, rdr1alleles, rdr2alleles, drloc)
    dqfprob = getAAgenostringmatchFIBERS(ddq1alleles,ddq2alleles,rdq1alleles,rdq2alleles, dqloc)

    dDRa1,dDRa2,rDRa1,rDRa2,DRcount,DRpos,DRpos1,DRend_pos = getAAgenostringmatch(dalldr1,dalldr2,ralldr1,ralldr2, drloc)
    dDQa1,dDQa2,rDQa1,rDQa2,DQcount,DQpos,DQpos1,DQend_pos = getAAgenostringmatch(dalldq1,dalldq2,ralldq1,ralldq2, dqloc)

    DRpos2, drprob = drfibersprob(DRpos)
    DQpos2, dqprob = dqfibersprob(DQpos)

    
    ddr1prob = dfdr1/ddr1sumant
    ddr2prob = dfdr2/ddr2sumant
    rdr1prob = rfdr1/rdr1sumant
    rdr2prob = rfdr2/rdr2sumant

    ddq1prob = dfdq1/ddq1sumant
    ddq2prob = dfdq2/ddq2sumant
    rdq1prob = rfdq1/rdq1sumant
    rdq2prob = rfdq2/rdq2sumant

    

    drrange = DRend_pos
    dqrange = DQend_pos

    #fprob = fibersprobability(locus, pos)
    prob = dqprob + drprob

    if (prob >0):
        fpred = "Yes"
        fhazard = 1.10
    else:
        fpred = "No"

    return render(request, 'DRQantigen2aa_out.html', 
        {'ddr1prob':ddr1prob, 'ddr2prob' : ddr2prob,'rdr1prob' : rdr1prob, 'rdr2prob' : rdr2prob, 'ddq1prob' : ddq1prob, 'ddq2prob' : ddq2prob, 'rdq1prob' : rdq1prob, 'rdq2prob' : rdq2prob, 'drace': drace,'dantdr1': dantdr1, 'dalldr1': dalldr1, 'dfdr1': dfdr1, 'dantdr2': dantdr2, 'dalldr2': dalldr2, 'dfdr2' : dfdr2, 'dantdq1': dantdq1, 'dalldq1':dalldq1,'dfdq1':dfdq1,'dantdq2':dantdq2,'dalldq2':dalldq2,'dfdq2':dfdq2, 'rrace': rrace, 'rantdr1': rantdr1, 'ralldr1': ralldr1, 'rfdr1': rfdr1, 'rantdr2': rantdr2, 'ralldr2': ralldr2, 'rfdr2' : rfdr2,'rantdq1': rantdq1, 'ralldq1':ralldq1,'rfdq1':rfdq1,'rantdq2':rantdq2,'ralldq2':ralldq2,'rfdq2':rfdq2, 'drrange': drrange, 'DRcount':DRcount, 'DRpos' : DRpos, 'DRpos1': DRpos1, 'DRpos2' : DRpos2, 'dqrange': dqrange, 'DQcount':DQcount, 'DQpos' : DQpos, 'DQpos1': DQpos1, 'DQpos2' : DQpos2, 'fpred' : fpred,'dqfprob' : dqfprob, 'drfprob': drfprob,'fhazard': fhazard})
