import numpy as np
import ROOT 
from ROOT import gStyle,gPad,RooFit,RooAbsData
import uproot 
import pandas as pd
from rootpy.vector import Vector2, Vector3, LorentzVector, Rotation, LorentzRotation
#import matplotlib.pyplot as plt

df=pd.read_table('Data_sig_tos_weights.txt',sep="\s+")

### generate plots and add data ###
def pltini(no,bin,data,xaxis,color,weight):
    h = ROOT.TH1F(no, "", bin, min(data), max(data))
    for m in range (len(data)):
	    h.Fill(data[m],weight[m])
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

def boost_to_mother(particle1,particle2,particle3,particle4,motherparicle):
    boosttomother = -(motherparicle.BoostVector())
    particle1.Boost(boosttomother)
    particle2.Boost(boosttomother)
    particle3.Boost(boosttomother)
    particle4.Boost(boosttomother)
    #return particle1, particle2, particle3, particle4


   
#========================== MAIN BODY
C_T_positive_list = []
C_T_negative_list = []
invmass_D0Dbar0KpPim_list_p = []
invmass_D0Dbar0KpPim_list_n = []
weight_list_p = []
weight_list_n = []

for i in range (len(df.index)):
    #extract database to form four-vectors (px,py,pz,E)
    locals()['D0_lv_%d'      % (i)] = LorentzVector(df.D0_PX[i],    df.D0_PY[i],    df.D0_PZ[i],    df.D0_PE[i])
    locals()['Dbar0_lv_%d'   % (i)] = LorentzVector(df.D0bar_PX[i], df.D0bar_PY[i], df.D0bar_PZ[i], df.D0bar_PE[i])
    locals()['Kplus_lv_%d'   % (i)] = LorentzVector(df.K_Kst0_PX[i],    df.K_Kst0_PY[i],    df.K_Kst0_PZ[i],    df.K_Kst0_PE[i])
    locals()['Piminus_lv_%d' % (i)] = LorentzVector(df.Pi_Kst0_PX[i],   df.Pi_Kst0_PY[i],   df.Pi_Kst0_PZ[i],   df.Pi_Kst0_PE[i])
    locals()['B0_lv_%d'      % (i)] = LorentzVector(df.B0_PX[i],   df.B0_PY[i],   df.B0_PZ[i],   df.B0_PE[i])
    locals()['sWeights_%d'   % (i)] = df.sWeights[i]
    
    boost_to_mother(locals()['D0_lv_%d'% (i)],locals()['Dbar0_lv_%d' % (i)],locals()['Kplus_lv_%d' % (i)],locals()['Piminus_lv_%d' % (i)],locals()['B0_lv_%d' % (i)])
    
    #check if it boost correctly
    locals()['B0_lv_%d' % (i)].Boost(-locals()['B0_lv_%d' % (i)].BoostVector())
    #print(locals()['B0_lv_%d' % (i)])
    
    #the sum of the four-momentum of the daughter-pair
    locals()['momD0Dbar0KpPim_%d' % (i)] = locals()['D0_lv_%d' % (i)] + locals()['Dbar0_lv_%d' % (i)] + locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
       
    #find the triple products for K+ (D0 X Dbar0)
    locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)])
    
    #find the invariant mass for the four particles
    locals()['invmass_D0Dbar0KpPim_%d' % (i)] = locals()['momD0Dbar0KpPim_%d' % (i)].M()


    if locals()['C_Kpd_D0cDbar0_%d' % (i)] > 0:
        C_T_positive_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])
        invmass_D0Dbar0KpPim_list_p.append(locals()['invmass_D0Dbar0KpPim_%d' % (i)])
        weight_list_p.append(locals()['sWeights_%d' % (i)])
    else:
        C_T_negative_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])
        invmass_D0Dbar0KpPim_list_n.append(locals()['invmass_D0Dbar0KpPim_%d' % (i)])
        weight_list_n.append(locals()['sWeights_%d' % (i)])



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

#Set up model
x_fit = ROOT.RooRealVar("x_fit","x_fit",min(invmass_D0Dbar0KpPim_list_p),max(invmass_D0Dbar0KpPim_list_p))
x = ROOT.RooRealVar("x","x",min(invmass_D0Dbar0KpPim_list_p),max(invmass_D0Dbar0KpPim_list_p))
mean = ROOT.RooRealVar("mean","mean of the gaussian",5200,5100,5300)
sigma = ROOT.RooRealVar("sigma","width of the gaussian",10,0.001,1000)
gauss = ROOT.RooGaussian("gauss","gaussian PDF",x_fit,mean,sigma)
w = ROOT.RooRealVar("w","w",min(weight_list_p),max(weight_list_p))

#varx = ROOT.RooRealVar("varx","varxdata",min(invmass_D0Dbar0KpPim_list_p),max(invmass_D0Dbar0KpPim_list_p))
data = ROOT.RooDataSet("data","dataset",ROOT.RooArgSet(x,w))
#data_w = ROOT.RooDataSet("data_w","weight data set",ROOT.RooArgSet(w))

for i in range (len(invmass_D0Dbar0KpPim_list_p)):
    x.setVal(invmass_D0Dbar0KpPim_list_p[i])
    w.setVal(weight_list_p[i])
    #w.Print("")
    data.add(ROOT.RooArgSet(x,w))
    #data_w.add(ROOT.RooArgSet(w))

wdata = ROOT.RooDataSet("data","dataset",data,ROOT.RooArgSet(x,w),"","w")
# varx = ROOT.RooRealVar("varx","varxdata",min(invmass_D0Dbar0KpPim_list_p),max(invmass_D0Dbar0KpPim_list_p))
# data = ROOT.RooDataHist("data","dataset",ROOT.RooArgSet(varx))
# for i in range (len(invmass_D0Dbar0KpPim_list_p)):
#     varx.setVal(invmass_D0Dbar0KpPim_list_p[i])
#     data.add(ROOT.RooArgSet(varx))
#data.Print("t")  
#data2 = ROOT.RooDataHist("data2","dataset2",ROOT.RooArgSet(data))

#xframe1 = varx.frame(ROOT.RooFit.Name("xframe"),ROOT.RooFit.Title("Gaussian p.d.f"),ROOT.RooFit.Bins(25))
#xframe = x.frame(ROOT.RooFit.Name("xframe"),ROOT.RooFit.Title("Gaussian p.d.f"),ROOT.RooFit.Bins(25))

#data.plotOn(xframe,ROOT.RooFit.DrawOption("B"))
#gauss.plotOn(xframe)
NS = ROOT.RooRealVar("NS","NS",50,1,10000)
x_fit.setRange("GaussRange",5100,5300)
gaussE  = ROOT.RooExtendPdf("gaussE","gaussE",gauss,NS,"GaussRange")
gaussE.fitTo(wdata,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(ROOT.kTRUE))

#rootfit_result = gauss.fitTo(wdata,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(ROOT.kTRUE))
xframe = x_fit.frame(ROOT.RooFit.Name("xframe"),ROOT.RooFit.Title("Gaussian p.d.f"))
gaussE.plotOn(xframe)
wdata.plotOn(xframe,RooFit.DataError(RooAbsData.SumW2))

mean.Print()
sigma.Print()
c = ROOT.TCanvas("c6", "root", 120, 120, 510, 500) 
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.6)
#xframe1.Draw()
xframe.Draw()
"""
c6 = ROOT.TCanvas("c6", "root", 120, 120, 510, 500) 
h6_1 = pltini("h6_1",bin,invmass_D0Dbar0KpPim_list_p,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kAzure+2,weight_list_p)
h6_2 = pltini("h6_2",bin,invmass_D0Dbar0KpPim_list_n,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kRed,weight_list_n)
leg6_1 = ROOT.TLegend(0.5, 0.82, 0.75, 0.89); leg6_1.AddEntry(h6_1,"\ B^0(C_T > 0)", "lep"); legset(leg6_1)
leg6_2 = ROOT.TLegend(0.5, 0.82, 0.75, 0.89); leg6_2.AddEntry(h6_2,"\ B^0(C_T < 0)", "lep"); legset(leg6_2)
h6_2.SetMaximum(1.3*h6_1.GetMaximum()); h6_1.SetMaximum(1.3*h6_1.GetMaximum())
pad6_1 = ROOT.TPad("pad6_1", "", 0, 0, 0.5, 1); padset(pad6_1); h6_1.Draw(); leg6_1.Draw();c6.cd()
pad6_2 = ROOT.TPad("pad6_2", "", 0.5, 0, 1, 1); padset(pad6_2); h6_2.Draw(); leg6_2.Draw()

#c6.SaveAs("results/B0_invmass_D0Dbar0KpPim_10000.png")
"""
raw_input()

