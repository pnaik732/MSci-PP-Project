import numpy as np
import ROOT 
from ROOT import gStyle,gPad
import uproot 
import pandas as pd
from rootpy.vector import Vector2, Vector3, LorentzVector, Rotation, LorentzRotation
#import matplotlib.pyplot as plt

df=pd.read_table('B0toDDbar0K+pi-_conj.txt',sep="\s+")

class BreitWigner:
    def __call__(self,t,par):
        M = par[0]
        Gamma = par[1]
        A = par[2]
        E = t[0]
        gamma = (M**2*(M**2+Gamma**2))**0.5
        k = (2**1.5*M*Gamma*gamma)/(np.pi*(M**2+gamma)**0.5)
        f = A * k/((E**2-M**2)**2+(M*Gamma)**2)
        return f

class BreitWigner1:
    def __call__(self,t,par):
        M1 = par[0]; Gamma1 = par[1]; A1 = par[2];
        M2 = par[3]; Gamma2 = par[4]; A2 = par[5];
        M3 = par[6]; Gamma3 = par[7]; A3 = par[8];
        E = t[0]
        f1 = BW0 (M1,Gamma1,A1,E)
        f2 = BW0 (M2,Gamma2,A2,E)
        f3 = BW0 (M3,Gamma3,A3,E)
        f = f1+f2+f3
        return f
    
def BW0 (M,Gamma,A,E):
    gamma = (M**2*(M**2+Gamma**2))**0.5
    k = (2**1.5*M*Gamma*gamma)/(np.pi*(M**2+gamma)**0.5)
    f = A * k/((E**2-M**2)**2+(M*Gamma)**2)
    return f

### generate plots and add data ###
def pltini(no,bin,data,xaxis,color):
    h = ROOT.TH1F(no, "", bin, min(data), max(data))
    for m in range (len(data)):
	    h.Fill(data[m])
    pltset (h,xaxis,color)
    return h

### customise the front size and colour ###    
def pltset (input, xaxis,color):
  #function to set up the fontsize of the axis labels and titles.
  gr = input
  gr.SetStats(False)
  gr.SetLineColor(color)
  gr.GetXaxis().SetTitle(xaxis)
  gr.GetYaxis().SetTitle('Candidates')
  gr.GetXaxis().SetLabelFont(132)
  gr.GetYaxis().SetLabelFont(132)
  gr.GetXaxis().SetTitleFont(132)
  gr.GetYaxis().SetTitleFont(132)
  gr.GetXaxis().SetLabelSize(0.05)
  gr.GetXaxis().SetTitleSize(0.05)
  gr.GetYaxis().SetLabelSize(0.05)
  gr.GetYaxis().SetTitleSize(0.05)
  gr.SetMinimum(0)

### add legend to the plot ###
def legset(input):
    legend = input
    legend.SetTextFont(133)
    legend.SetTextSize(20)
    legend.SetLineColor(0)

def fit(left,right,color,M,Gamma):
    fit = ROOT.TF1("f",BreitWigner(),left,right,3)
    fit.SetLineColor(color)
    fit.SetParameters(M,Gamma,10)
    return fit

def fit1(left,right,M1,Gamma1,A1, M2,Gamma2,A2, M3,Gamma3,A3):
    fit = ROOT.TF1("f",BreitWigner1(),left,right,9)
    fit.SetParameters(M1,Gamma1,A1, M2,Gamma2,A2, M3,Gamma3,A3)
    return fit

def padset(input):
    pad = input
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.1)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.15)
    pad.Draw()
    pad.cd()

### find the triple product in the rest frame of mother particle###
def tripleproduct(particleA1,particleA2,particleB):
    particleA1v = particleA1.Vect()
    particleA2v = particleA2.Vect()
    particleBv  = particleB.Vect()
    cross       = particleA1v.Cross(particleA2v)
    dot         = particleBv.Dot(cross)
    return dot

### find the helicity angle in the rest frame  of daughter pairs###
def coshel(particle, parent, motherparticle):
    boosttoparent = -(parent.BoostVector())
    particle.Boost(boosttoparent)
    motherparticle.Boost(boosttoparent)
    particle3 = particle.Vect()
    motherparticle3 = motherparticle.Vect()
    numerator = particle3.Dot(motherparticle3)
    denominator = (particle3.Mag()) * (motherparticle3.Mag())
    cosTheta = numerator/denominator
    return cosTheta

### find the angle between the two planes in the rest frame of mother particle###
def planeangle(particleA1,particleA2,particleB1,particleB2,motherparticle):
    particleA1v = particleA1.Vect()
    particleA2v = particleA2.Vect()
    particleB1v = particleB1.Vect()
    particleB2v = particleB2.Vect()
    particleAv  = (particleA1 + particleA2).Vect()
    phi         = (particleA1v.Cross(particleA2v)).Angle(particleB1v.Cross(particleB2v))
    #cross = ((particleA1v.Cross(particleA2v)) .Cross((particleB1v.Cross(particleB2v)) )).Angle(particleAv)
    #if  cross < np.pi*0.5:
    #    phi = - phi
    return phi
   
#========================== MAIN BODY
invmass_D0Dbar0_list = []
invmass_D0Dbar0_list_p = []
invmass_D0Dbar0_list_n = []
invmass_KpPim_list = []
invmass_KpPim_list_p = []
invmass_KpPim_list_n = []
cosTheta_D0Dbar0_D0_list_p = []
cosTheta_D0Dbar0_D0_list_n = []
cosTheta_KpPim_Kp_list_p = []
cosTheta_KpPim_Kp_list_n = []
Phi_list_p = []
Phi_list_n = []
C_T_positive_list = []
C_T_negative_list = []
#invmass_D0Dbar0KpPim_list_p = []
#invmass_D0Dbar0KpPim_list_n = []

for i in range (len(df.index)):
    #extract database to form four-vectors (px,py,pz,E)
    locals()['D0_lv_%d'      % (i)] = LorentzVector(-df._1_D0_Px[i],    -df._1_D0_Py[i],    -df._1_D0_Pz[i],    -df._1_D0_E[i])
    locals()['Dbar0_lv_%d'   % (i)] = LorentzVector(-df._2_Dbar0_Px[i], -df._2_Dbar0_Py[i], -df._2_Dbar0_Pz[i], -df._2_Dbar0_E[i])
    locals()['Kplus_lv_%d'   % (i)] = LorentzVector(-df._3_Kp_Px[i],    -df._3_Kp_Py[i],    -df._3_Kp_Pz[i],    -df._3_Kp_E[i])
    locals()['Piminus_lv_%d' % (i)] = LorentzVector(-df._4_pim_Px[i],   -df._4_pim_Py[i],   -df._4_pim_Pz[i],   -df._4_pim_E[i])
    
    #the sum of the four-momentum of the daughter-pair
    locals()['momD0Dbar0_%d' % (i)] = locals()['D0_lv_%d'    % (i)] + locals()['Dbar0_lv_%d'   % (i)]
    locals()['momKpPim_%d'   % (i)] = locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
    locals()['momD0Dbar0KpPim_%d' % (i)] = locals()['D0_lv_%d' % (i)] + locals()['Dbar0_lv_%d' % (i)] + locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]

    #find the five CM variables
    locals()['invmass_D0Dbar0_%d'     % (i)] = locals()['momD0Dbar0_%d' % (i)].M()
    locals()['invmass_KpPim_%d'       % (i)] = locals()['momKpPim_%d'   % (i)].M()
    locals()['cosTheta_D0Dbar0_D0_%d' % (i)] = coshel(locals()['D0_lv_%d'    % (i)], locals()['momD0Dbar0_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
    locals()['cosTheta_KpPim_Kp_%d'   % (i)] = coshel(locals()['Kplus_lv_%d' % (i)], locals()['momKpPim_%d'   % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
    locals()['Phi_%d'                 % (i)] = planeangle(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)], locals()['Piminus_lv_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
    
    #find the triple products for K+ (D0 X Dbar0)
    locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)])
    
    #find the invariant mass for the four particles
    #locals()['invmass_D0Dbar0KpPim_%d' % (i)] = locals()['momD0Dbar0KpPim_%d' % (i)].M()

    #check the if the boost is correct
    #print((-locals()['momD0Dbar0_%d' % (i)].BoostVector()).Y())
    
    #locals()['momD0Dbar0_%d' % (i)].Boost(-locals()['momD0Dbar0_%d' % (i)].BoostVector())
    #print(locals()['momD0Dbar0_%d' % (i)])
    invmass_D0Dbar0_list.append(locals()['invmass_D0Dbar0_%d' % (i)])
    invmass_KpPim_list.append(locals()['invmass_KpPim_%d'   % (i)])
    if -locals()['C_Kpd_D0cDbar0_%d' % (i)] > 0:
        invmass_D0Dbar0_list_p.append(locals()['invmass_D0Dbar0_%d' % (i)])
        invmass_KpPim_list_p.append(locals()['invmass_KpPim_%d'   % (i)])
        cosTheta_D0Dbar0_D0_list_p.append(locals()['cosTheta_D0Dbar0_D0_%d' % (i)])
        cosTheta_KpPim_Kp_list_p.append(locals()['cosTheta_KpPim_Kp_%d'   % (i)])
        Phi_list_p.append(locals()['Phi_%d' % (i)])
        C_T_positive_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])
        #invmass_D0Dbar0KpPim_list_p.append(locals()['invmass_D0Dbar0KpPim_%d' % (i)])
    else:
        invmass_D0Dbar0_list_n.append(locals()['invmass_D0Dbar0_%d' % (i)])
        invmass_KpPim_list_n.append(locals()['invmass_KpPim_%d' % (i)])
        cosTheta_D0Dbar0_D0_list_n.append(locals()['cosTheta_D0Dbar0_D0_%d' % (i)])
        cosTheta_KpPim_Kp_list_n.append(locals()['cosTheta_KpPim_Kp_%d' % (i)])
        Phi_list_n.append(locals()['Phi_%d' % (i)])
        C_T_negative_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])
        #invmass_D0Dbar0KpPim_list_n.append(locals()['invmass_D0Dbar0KpPim_%d' % (i)])

#Triple product asymmetry and the error
C_T_p = 1.0 * len(C_T_positive_list)
C_T_n = 1.0 * len(C_T_negative_list)
Asy   = (C_T_p-C_T_n) / (C_T_p+C_T_n)
err_A = ((2 * C_T_n / (C_T_p + C_T_n) ** 2) ** 2 * C_T_p + (2 * C_T_p / (C_T_p + C_T_n) ** 2) ** 2 * C_T_n) ** 0.5
print(Asy, err_A)

#draw the distributions of CM variables
#gStyle.SetPadLeftMargin(0.12)
#gStyle.SetPadBottomMargin(0.12);
bin = 100
A1,A2,A3 = 10,10,10

c0 = ROOT.TCanvas("c0", "root", 120, 120, 1100, 500) 
h0_1 = pltini("h0_1",bin,invmass_D0Dbar0_list,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kBlack)
h0_2 = pltini("h0_2",bin,invmass_KpPim_list,  '\ m(K^+\pi^-)\\;[GeV/c^2]',     ROOT.kBlack)
f0_1 = fit1(min(invmass_D0Dbar0_list),max(invmass_D0Dbar0_list),3.770,0.0272,A1,4.039,0.080,A2,4.191,0.070,A3); h0_1.Fit(f0_1,"R+")
f0_2 = fit(0.700,1.400,4,0.89555,0.0473); h0_2.Fit(f0_2,"R+")
leg0_1 = ROOT.TLegend(0.66, 0.65, 0.89, 0.89); leg0_1.AddEntry(h0_1,"\ m(D^0\\bar{D}^0)", "lep"); leg0_1.AddEntry(f0_1,"fit", "l"); legset(leg0_1)
leg0_2 = ROOT.TLegend(0.66, 0.65, 0.89, 0.89); leg0_2.AddEntry(h0_2,"\ m(K^+\pi^-)", "lep");      leg0_2.AddEntry(f0_2,"fit", "l"); legset(leg0_2)
pad0_1 = ROOT.TPad("pad0_1", "", 0, 0, 0.5, 1); padset(pad0_1); h0_1.Draw(); leg0_1.Draw();c0.cd()
pad0_2 = ROOT.TPad("pad0_2", "", 0.5, 0, 1, 1); padset(pad0_2); h0_2.Draw(); leg0_2.Draw()

c0.SaveAs("results/B0_conj2_invariant_mass_fit_10000.png")

c1 = ROOT.TCanvas("c1", "root", 120, 120, 1100, 500) 
h1_1 = pltini("h1_1",bin,invmass_D0Dbar0_list_p,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kAzure+2)
h1_2 = pltini("h1_2",bin,invmass_D0Dbar0_list_n,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kRed)
f1_1 = fit1(min(invmass_D0Dbar0_list_p),max(invmass_D0Dbar0_list_p),3.770,0.0272,A1,4.039,0.080,A2,4.191,0.070,A3);h1_1.Fit(f1_1,"R+")
f1_2 = fit1(min(invmass_D0Dbar0_list_n),max(invmass_D0Dbar0_list_n),3.770,0.0272,A1,4.039,0.080,A2,4.191,0.070,A3);h1_2.Fit(f1_2,"R+")
leg1_1 = ROOT.TLegend(0.46, 0.65, 0.75, 0.89); leg1_1.AddEntry(h1_1,"\ B^0(C_T > 0)", "lep"); leg1_1.AddEntry(f1_1,"fit", "l"); legset(leg1_1)
leg1_2 = ROOT.TLegend(0.46, 0.65, 0.75, 0.89); leg1_2.AddEntry(h1_2,"\ B^0(C_T < 0)", "lep"); leg1_2.AddEntry(f1_2,"fit", "l"); legset(leg1_2)
h1_2.SetMaximum(1.3*h1_1.GetMaximum()); h1_1.SetMaximum(1.3*h1_1.GetMaximum())
pad1_1 = ROOT.TPad("pad1_1", "", 0, 0, 0.5, 1); padset(pad1_1); h1_1.Draw(); leg1_1.Draw();c1.cd()
pad1_2 = ROOT.TPad("pad1_2", "", 0.5, 0, 1, 1); padset(pad1_2); h1_2.Draw(); leg1_2.Draw()

c1.SaveAs("results/B0_conj2_invmass_D0Dbar0_fit_10000.png")

c2 = ROOT.TCanvas("c2", "root", 120, 120, 1100, 500) 
h2_1 = pltini("h2_1",bin,invmass_KpPim_list_p,'\ m(K^+\pi^-)\\;[GeV/c^2]',ROOT.kAzure+2)
h2_2 = pltini("h2_2",bin,invmass_KpPim_list_n,'\ m(K^+\pi^-)\\;[GeV/c^2]',ROOT.kRed)
f2_1 = fit(0.700,1.400,4,0.89176,0.0473); h2_1.Fit(f2_1,"R+")
f2_2 = fit(0.700,1.400,4,0.89176,0.0473); h2_2.Fit(f2_2,"R+")
leg2_1 = ROOT.TLegend(0.56, 0.8, 0.75, 0.89); leg2_1.AddEntry(h2_1,"\ B^0(C_T > 0)", "lep"); leg2_1.AddEntry(f2_1,"\ K^{*0}(892)", "l"); legset(leg2_1)
leg2_2 = ROOT.TLegend(0.56, 0.8, 0.75, 0.89); leg2_2.AddEntry(h2_2,"\ B^0(C_T < 0)", "lep"); leg2_2.AddEntry(f2_2,"\ K^{*0}(892)", "l"); legset(leg2_2)
h2_2.SetMaximum(1.3*h2_1.GetMaximum()); h2_1.SetMaximum(1.3*h2_1.GetMaximum())
pad2_1 = ROOT.TPad("pad2_1", "", 0, 0, 0.5, 1); padset(pad2_1); h2_1.Draw(); leg2_1.Draw();c2.cd()
pad2_2 = ROOT.TPad("pad2_2", "", 0.5, 0, 1, 1); padset(pad2_2); h2_2.Draw(); leg2_2.Draw()
    
c2.SaveAs("results/B0_conj2_invmass_KpPim_fit_10000.png")

with open('results/data_conj.txt',mode='a') as f:
    f.write("bin=%d A1=%d A2=%d A3=%d"%(bin,A1,A2,A3)+ "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(3770)",f0_1.GetParameters()[0],f0_1.GetParError(0), f0_1.GetParameters()[1],f0_1.GetParError(1),f0_1.GetChisquare(),"combined") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(3770)",f1_1.GetParameters()[0],f1_1.GetParError(0), f1_1.GetParameters()[1],f1_1.GetParError(1),f1_1.GetChisquare(),"C_T > 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(3770)",f1_2.GetParameters()[0],f1_2.GetParError(0), f1_2.GetParameters()[1],f1_2.GetParError(1),f1_2.GetChisquare(),"C_T < 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4040)",f0_1.GetParameters()[3],f0_1.GetParError(3), f0_1.GetParameters()[4],f0_1.GetParError(4),f0_1.GetChisquare(),"combined") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4040)",f1_1.GetParameters()[3],f1_1.GetParError(3), f1_1.GetParameters()[4],f1_1.GetParError(4),f1_1.GetChisquare(),"C_T > 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4040)",f1_2.GetParameters()[3],f1_2.GetParError(3), f1_2.GetParameters()[4],f1_2.GetParError(4),f1_2.GetChisquare(),"C_T < 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4160)",f0_1.GetParameters()[6],f0_1.GetParError(6), f0_1.GetParameters()[7],f0_1.GetParError(7),f0_1.GetChisquare(),"combined") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4160)",f1_1.GetParameters()[6],f1_1.GetParError(6), f1_1.GetParameters()[7],f1_1.GetParError(7),f1_1.GetChisquare(),"C_T > 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4160)",f1_2.GetParameters()[6],f1_2.GetParError(6), f1_2.GetParameters()[7],f1_2.GetParError(7),f1_2.GetChisquare(),"C_T < 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("K*0(892)", f0_2.GetParameters()[0],f0_2.GetParError(0), f0_2.GetParameters()[1],f0_2.GetParError(1),f0_2.GetChisquare(),"combined") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("K*0(892)", f2_1.GetParameters()[0],f2_1.GetParError(0), f2_1.GetParameters()[1],f2_1.GetParError(1),f2_1.GetChisquare(),"C_T > 0") + "\n")
    f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("K*0(892)", f2_2.GetParameters()[0],f2_2.GetParError(0), f2_2.GetParameters()[1],f2_2.GetParError(1),f2_2.GetChisquare(),"C_T < 0") + "\n")


c3 = ROOT.TCanvas("c3", "root", 120, 120, 1100, 500) 
h3_1 = pltini("h3_1",bin,cosTheta_D0Dbar0_D0_list_p,'\ \cos(\\theta_D)',ROOT.kAzure+2)
h3_2 = pltini("h3_2",bin,cosTheta_D0Dbar0_D0_list_n,'\ \cos(\\theta_D)',ROOT.kRed)
leg3_1 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg3_1.AddEntry(h3_1,"\ B^0(C_T > 0)", "lep"); legset(leg3_1)
leg3_2 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg3_2.AddEntry(h3_2,"\ B^0(C_T < 0)", "lep"); legset(leg3_2)
pad3_1 = ROOT.TPad("pad3_1", "", 0, 0, 0.5, 1); padset(pad3_1); h3_1.Draw(); leg3_1.Draw();c3.cd()
h3_2.SetMaximum(1.3*h3_1.GetMaximum()); h3_1.SetMaximum(1.3*h3_1.GetMaximum())
pad3_2 = ROOT.TPad("pad3_2", "", 0.5, 0, 1, 1); padset(pad3_2); h3_2.Draw(); leg3_2.Draw()

c3.SaveAs("results/B0_conj2_helangle_D0Dbar0_10000.png")

c4 = ROOT.TCanvas("c4", "root", 120, 120, 1100, 500) 
h4_1 = pltini("h4_1",bin,cosTheta_KpPim_Kp_list_p,'\ \cos(\\theta_K)',ROOT.kAzure+2)
h4_2 = pltini("h4_2",bin,cosTheta_KpPim_Kp_list_n,'\ \cos(\\theta_K)',ROOT.kRed)
leg4_1 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg4_1.AddEntry(h4_1,"\ B^0(C_T > 0)", "lep"); legset(leg4_1)
leg4_2 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg4_2.AddEntry(h4_2,"\ B^0(C_T < 0)", "lep"); legset(leg4_2)
h4_2.SetMaximum(1.3*h4_1.GetMaximum()); h4_1.SetMaximum(1.3*h4_1.GetMaximum())
pad4_1 = ROOT.TPad("pad4_1", "", 0, 0, 0.5, 1); padset(pad4_1); h4_1.Draw(); leg4_1.Draw();c4.cd()
pad4_2 = ROOT.TPad("pad4_2", "", 0.5, 0, 1, 1); padset(pad4_2); h4_2.Draw(); leg4_2.Draw()

c4.SaveAs("results/B0_conj2_helangle_KpPim_10000.png")

c5 = ROOT.TCanvas("c5", "root", 120, 120, 1100, 500) 
h5_1 = pltini("h5_1",bin,Phi_list_p,'\ \phi\\;[rad]',ROOT.kAzure+2)
h5_2 = pltini("h5_2",bin,Phi_list_n,'\ \phi\\;[rad]',ROOT.kRed)
leg5_1 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg5_1.AddEntry(h5_1,"\ B^0(C_T > 0)", "lep"); legset(leg5_1)
leg5_2 = ROOT.TLegend(0.56, 0.82, 0.75, 0.89); leg5_2.AddEntry(h5_2,"\ B^0(C_T < 0)", "lep"); legset(leg5_2)
h5_2.SetMaximum(1.3*h5_1.GetMaximum()); h5_1.SetMaximum(1.3*h5_1.GetMaximum())
pad5_1 = ROOT.TPad("pad5_1", "", 0, 0, 0.5, 1); padset(pad5_1); h5_1.Draw(); leg5_1.Draw();c5.cd()
pad5_2 = ROOT.TPad("pad5_2", "", 0.5, 0, 1, 1); padset(pad5_2); h5_2.Draw(); leg5_2.Draw()

c5.SaveAs("results/B0_conj2_planeangle_10000.png")

#c6 = ROOT.TCanvas("c6", "root", 120, 120, 510, 500) 
#h6_1 = pltini("h6_1",bin,invmass_D0Dbar0KpPim_list_p,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kAzure+2)
#h6_2 = pltini("h6_2",bin,invmass_D0Dbar0KpPim_list_n,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kRed)
#leg5_1 = ROOT.TLegend(0.5, 0.82, 0.75, 0.89); leg5_1.AddEntry(h5_1,"\ B^0(C_T > 0)", "lep"); legset(leg5_1)
#leg5_2 = ROOT.TLegend(0.5, 0.82, 0.75, 0.89); leg5_2.AddEntry(h5_2,"\ B^0(C_T < 0)", "lep"); legset(leg5_2)
#h6_2.SetMaximum(1.3*h6_1.GetMaximum()); h6_1.SetMaximum(1.3*h6_1.GetMaximum())
#pad5_1 = ROOT.TPad("pad5_1", "", 0, 0, 0.5, 1); padset(pad5_1); h5_1.Draw(); leg5_1.Draw();c5.cd()
#pad5_2 = ROOT.TPad("pad5_2", "", 0.5, 0, 1, 1); padset(pad5_2); h5_2.Draw(); leg5_2.Draw()

#c6.SaveAs("results/B0_conj2_invmass_D0Dbar0KpPim_10000.png")

raw_input()

