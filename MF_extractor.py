import pylhe
import hist
import numpy as np

b_ind=-1
nu_ind=-2
tau_ind=-3

def a_histo(nbs,vmin,vmax,values,xlab,ylab="Count",log_flag=False):
    HistA = hist.Hist.new.Reg(nbs, vmin, vmax).Int64()
    HistA.fill(values)
    HistAart = HistA.plot1d()
    HistAax = HistAart[0].stairs.axes
    if log_flag: HistAax.set_yscale("log")
    HistAax.set_xlabel(xlab)
    HistAax.set_ylabel(ylab)
    return HistA, HistAart

def pt_hist(myevts,nbins,ptmin,ptmax,part_id,part_name,mlog=False):
    return a_histo(nbins,ptmin,ptmax,
                   (myevts.particles.vector[:, part_id]).pt,
                   r"p_T("+part_name+") [GeV]",
                   log_flag=mlog)

def eta_hist(myevts,nbins,etamin,etamax,part_id,part_name,mlog=False):
    return a_histo(nbins, etamin, etamax,
                   (myevts.particles.vector[:, part_id]).eta,
                   r"\eta("+part_name+")",
                   log_flag=mlog)

def phi_hist(myevts,nbins,phimin,phimax,part_id,part_name,mlog=False):
    return a_histo(nbins, phimin, phimax,
                   (myevts.particles.vector[:, part_id]).phi,
                   r"\phi("+part_name+")",
                   log_flag=mlog)

def minv_hist(myevts,nbins,minvmin,minvmax,p1_id,p2_id,p1_name,p2_name,mlog=False):
    return a_histo(nbins, minvmin, minvmax,
                   (myevts.particles.vector[:, p1_id] + myevts.particles.vector[:, p2_id]).mass,
                   r"m("+p1_name+","+p2_name+") [GeV]",
                   log_flag=mlog)

def mT_nu_hist(myevts,nbins,mTmin,mTmax,p1_id,p2_id,p1_name,mlog=False):
    return a_histo(nbins, mTmin, mTmax,
                   np.sqrt((2*(myevts.particles.vector[:, p1_id]).pt*(myevts.particles.vector[:, p2_id]).pt)*(1-np.cos(myevts.particles.vector[:, p1_id].deltaphi(myevts.particles.vector[:, p2_id])+np.pi))),
                   r"mT("+p1_name+",p_{T}(\nu)) [GeV]",
                   log_flag=mlog)

def mTot_nu_hist(myevts,nbins,mTotmin,mTotmax,p1_id,p2_id,p3_id,p1_name,p2_name,mlog=False):
    FirstTerm=np.square(myevts.particles.vector[:, p1_id].pt+myevts.particles.vector[:, p2_id].pt+myevts.particles.vector[:, p3_id].pt)
    SecondTerm=np.square(((myevts.particles.vector[:, p1_id])+(myevts.particles.vector[:, p2_id])+(myevts.particles.vector[:, p3_id])).pt)
    return a_histo(nbins, mTotmin, mTotmax,
                   np.sqrt(FirstTerm-SecondTerm),
                   r"m_{tot}("+p1_name+","+p2_name+",p_{T}^{miss}) [GeV]",
                   log_flag=mlog)

def Phi_0_2pi(phiangle):
    mask1 = phiangle >= 2*np.pi
    mask2 = phiangle < 0
    newphiangle = phiangle+(mask1*-2*np.pi)+(mask2*2*np.pi)
    return newphiangle

def Met_Vec_2p(myevts,p1_id,p2_id):
    MET = np.sqrt( np.square(myevts.particles.vector[:, p1_id].px + myevts.particles.vector[:, p2_id].px) +np.square(myevts.particles.vector[:, p1_id].py + myevts.particles.vector[:, p2_id].py) )
    METphi=Phi_0_2pi(np.arctan( (myevts.particles.vector[:, p1_id].py + myevts.particles.vector[:, p2_id].py) / (myevts.particles.vector[:, p1_id].px + myevts.particles.vector[:, p2_id].px) )+2*np.pi)
    return MET, METphi

def mT_MET_hist(myevts,nbins,mTmin,mTmax,p1_id,p2_id,p1_name,mlog=False):
    MET, METphi = Met_Vec_2p(myevts,p1_id,p2_id)
    return a_histo(nbins, mTmin, mTmax,
                   np.sqrt((2*(myevts.particles.vector[:, p1_id]).pt*(MET)*(1-np.cos(Phi_0_2pi(Phi_0_2pi(myevts.particles.vector[:, p1_id].phi)-METphi))))),
                   r"mT("+p1_name+",p_{T}^{miss}) [GeV]",
                   log_flag=mlog)

def mTot_MET_hist(myevts,nbins,mTotmin,mTotmax,p1_id,p2_id,p1_name,p2_name,mlog=False):
    MET, METphi = Met_Vec_2p(myevts,p1_id,p2_id)
    FirstTerm=np.square(myevts.particles.vector[:, p1_id].pt+myevts.particles.vector[:, p2_id].pt+MET)
    SecondTerm=np.square(myevts.particles.vector[:, p1_id].px+myevts.particles.vector[:, p2_id].px+(MET*np.cos(METphi))) + np.square(myevts.particles.vector[:, p1_id].py+myevts.particles.vector[:, p2_id].py+(MET*np.sin(METphi))) 
    return a_histo(nbins, mTotmin, mTotmax,
                   np.sqrt(FirstTerm-SecondTerm),
                   r"m_{tot}("+p1_name+","+p2_name+",p_{T}^{miss}) [GeV]",
                   log_flag=mlog)


########################################
#OLD
########################################

#def pt_hist(nbins,ptmin,ptmax,part_id,part_name,mlog=False):
#    HistPT = hist.Hist.new.Reg(nbins, ptmin, ptmax).Int64()
#    HistPT.fill((events.particles.vector[:, part_id]).pt)
#    HistPTart = HistPT.plot1d()
#    HistPTax = HistPTart[0].stairs.axes
#    if log_flag: HistPTax.set_yscale("log")
#    HistPTax.set_xlabel(r"p_T("+part_name+") [GeV]")
#    HistPTax.set_ylabel("Count")
#def eta_hist(nbins,etamin,etamax,part_id,part_name,log_flag=False):
#    HistETA = hist.Hist.new.Reg(nbins, etamin, etamax).Int64()
#    HistETA.fill((events.particles.vector[:, part_id]).eta)
#    HistETAart = HistETA.plot1d()
#    HistETAax = HistETAart[0].stairs.axes
#    if log_flag: HistETAax.set_yscale("log")
#    HistETAax.set_xlabel(r"\eta("+part_name+")")
#    HistETAax.set_ylabel("Count")
#
#def phi_hist(nbins,phimin,phimax,part_id,part_name,log_flag=False):
#    HistPHI = hist.Hist.new.Reg(nbins, phimin, phimax).Int64()
#    HistPHI.fill((events.particles.vector[:, part_id]).pt)
#    HistPHIart = HistPHI.plot1d()
#    HistPHIax = HistPHIart[0].stairs.axes
#    if log_flag: HistPHIax.set_yscale("log")
#    HistPHIax.set_xlabel("\phi("+part_name+")")
#    HistPHIax.set_ylabel("Count")
#
#def minv_hist(nbins,minvmin,minvmax,p1_id,p2_id,p1_name,p2_name,log_flag=False):
#    HistM = hist.Hist.new.Reg(nbins, ptmin, ptmax).Int64()
#    HistM.fill((events.particles.vector[:, part_id]).pt)
#    HistMart = HistM.plot1d()
#    HistMax = HistMart[0].stairs.axes
#    if log_flag: HistMax.set_yscale("log")
#    HistMax.set_xlabel("m("+p1_name+","+p2_name+") [GeV]")
#    HistMax.set_ylabel("Count")
