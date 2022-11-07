import sys
from tracemalloc import start
import requests
from django.shortcuts import render
from django.http import HttpResponse
from .aa_matching import getAAposition
from .aa_matching import isPositionMismatched
from .aa_fibers import antigen2HFallele,getAAgenostringmatch, drfibersprob, dqfibersprob


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

    dalldr1,dfdr1,dnalldr1, ddr1sumant,ddr1allelefreqs = antigen2HFallele(drace,dantdr1)
    dalldr2,dfdr2,dnalldr2, ddr2sumant,ddr2allelefreqs = antigen2HFallele(drace,dantdr2)

    dalldq1,dfdq1,dnalldq1, ddq1sumant, ddq1allelefreqs = antigen2HFallele(drace,dantdq1)
    dalldq2,dfdq2,dnalldq2, ddq2sumant, ddq2allelefreqs = antigen2HFallele(drace,dantdq2)


    #get Recipient most prob alleles

    ralldr1,rfdr1,rnalldr1, rdr1sumant,rdr1allelefreqs = antigen2HFallele(rrace,rantdr1)
    ralldr2,rfdr2,rnalldr2, rdr2sumant,rdr2allelefreqs = antigen2HFallele(rrace,rantdr2)

    ralldq1,rfdq1,rnalldq1, rdq1sumant,rdq1allelefreqs = antigen2HFallele(rrace,rantdq1)
    ralldq2,rfdq2,rnalldq2, rdq2sumant,rdq2allelefreqs= antigen2HFallele(rrace,rantdq2)


    dDRa1,dDRa2,rDRa1,rDRa2,DRcount,DRpos,DRpos1,DRend_pos = getAAgenostringmatch(dalldr1,dalldr2,ralldr1,ralldr2, drloc)
    dDQa1,dDQa2,rDQa1,rDQa2,DQcount,DQpos,DQpos1,DQend_pos = getAAgenostringmatch(dalldq1,dalldq2,ralldq1,ralldq2, dqloc)


    drrange = DRend_pos
    dqrange = DQend_pos


    DRpos2, drprob = drfibersprob(DRpos)
    DQpos2, dqprob = dqfibersprob(DQpos)
    #fprob = fibersprobability(locus, pos)
    prob = dqprob + drprob

    if (prob >0):
        fpred = "Yes"
        fhazard = 1.10

        if (drprob > 0):
            ddr1 = dfdr1/ddr1sumant
            ddr2 = dfdr2/ddr2sumant
            ddr = ddr1 + ddr2
            ddrprob = ddr/2


            rdr1 = rfdr1/rdr1sumant
            rdr2 = rfdr2/rdr2sumant
            rdr = rdr1 + rdr2
            rdrprob = rdr/2
            print("prob of R:", rdrprob)
            print("prob of D:", ddrprob)

            drprob = ddrprob + rdrprob
            drfprob = drprob/2
            print("prob of DR:", drfprob)



        else:
            dqfprob = 0

        if (dqprob > 0):
            ddq1 = dfdq1/ddq1sumant
            ddq2 = dfdq2/ddq2sumant
            ddq = ddq1 + ddq2
            ddqprob = ddq/2


            rdq1 = rfdq1/rdq1sumant
            rdq2 = rfdq2/rdq2sumant
            rdq = rdq1 + rdq2
            rdqprob = rdq/2
            print("prob of R:", rdqprob)
            print("prob of D:", ddqprob)

            dqprob = ddqprob + rdqprob
            dqfprob = dqprob/2
            print("prob of DR:", dqfprob)
            

        else:
            dqfprob = 0
    else:
        fpred = "No"

    return render(request, 'DRQantigen2aa_out.html', 
        {'drace': drace,'dantdr1': dantdr1, 'dalldr1': dalldr1, 'dfdr1': dfdr1, 'dantdr2': dantdr2, 'dalldr2': dalldr2, 'dfdr2' : dfdr2, 'dantdq1': dantdq1, 'dalldq1':dalldq1,'dfdq1':dfdq1,'dantdq2':dantdq2,'dalldq2':dalldq2,'dfdq2':dfdq2, 'rrace': rrace, 'rantdr1': rantdr1, 'ralldr1': ralldr1, 'rfdr1': rfdr1, 'rantdr2': rantdr2, 'ralldr2': ralldr2, 'rfdr2' : rfdr2,'rantdq1': rantdq1, 'ralldq1':ralldq1,'rfdq1':rfdq1,'rantdq2':rantdq2,'ralldq2':ralldq2,'rfdq2':rfdq2, 'drrange': drrange, 'DRcount':DRcount, 'DRpos' : DRpos, 'DRpos1': DRpos1, 'DRpos2' : DRpos2, 'dqrange': dqrange, 'DQcount':DQcount, 'DQpos' : DQpos, 'DQpos1': DQpos1, 'DQpos2' : DQpos2, 'fpred' : fpred,'dqfprob' : dqfprob, 'drfprob': drfprob,'fhazard': fhazard})
