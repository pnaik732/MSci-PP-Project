import numpy as np
import ROOT
from ROOT import gPad
import sys

class BreitWigner:
    def __call__(self,t,par):
        M = par[0]
        Gamma = par[1]
        A = par[2]
        E = t[0]
        gamma = (M**2*(M**2+Gamma**2))**0.5
        k = (2**1.8*M*Gamma*gamma)/(np.pi*(M**2+gamma)**0.5)
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
    k = (2**1.8*M*Gamma*gamma)/(np.pi*(M**2+gamma)**0.5)
    f = A * k/((E**2-M**2)**2+(M*Gamma)**2)
    return f

def fit(left,right,color,M,Gamma):
    fit = ROOT.TF1("f",BreitWigner(),left,right,3)
    fit.SetLineColor(color)
    fit.SetParameters(M,Gamma,10)
    return fit

def fit1(left,right,M1,Gamma1,A1, M2,Gamma2,A2, M3,Gamma3,A3):
    fit = ROOT.TF1("f",BreitWigner1(),left,right,9)
    fit.SetParameters(M1,Gamma1,A1, M2,Gamma2,A2, M3,Gamma3,A3)
    return fit
    
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
    legend.SetTextSize(25)
    legend.SetLineColor(0)

def main(program,type,event,phase,event_n):   
	CM_data      = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n)))
	CM_data_conj = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d_conj.txt'%(phase,event_n,type,event,event_n)))

	mDDbar      = CM_data[1];            mDDbar_conj      = CM_data_conj[1]
	mKpi        = CM_data[2];            mKpi_conj        = CM_data_conj[2]

	bin = 400
	A1,A2,A3 = 10,10,10

	c0 = ROOT.TCanvas("c0", "B0 decay invmass", 1420,710); c0.Divide(2,1) 
	c0.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h0_1 = pltini("h0_1",bin,mDDbar,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kBlack)
	f0_1 = fit1(min(mDDbar),max(mDDbar),3.770,0.0272,A1,4.039,0.080,A2,4.191,0.070,A3); h0_1.Fit(f0_1,"R+")
	leg0_1 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg0_1.AddEntry(h0_1,"\ B^0 \\; data \\;(MC)", "lf"); leg0_1.AddEntry(f0_1,"\ BW \\; fit", "l"); legset(leg0_1)
	h0_1.SetMaximum(1.8*h0_1.GetMaximum())
	h0_1.Draw(); leg0_1.Draw() 

	c0.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h0_2 = pltini("h0_2",bin,mDDbar_conj,'\ m(\\bar{D}^0 D^0)\\;[GeV/c^2]',ROOT.kBlack)
	f0_2 = fit1(min(mDDbar_conj),max(mDDbar_conj),3.770,0.0272,A1,4.039,0.080,A2,4.191,0.070,A3); h0_2.Fit(f0_2,"R+")
	leg0_2 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg0_2.AddEntry(h0_2,"\  \\bar{B}^0  \\; data \\;(MC)", "lf"); leg0_2.AddEntry(f0_2,"\ BW \\; fit", "l"); legset(leg0_2)
	h0_2.SetMaximum(h0_1.GetMaximum())
	h0_2.Draw(); leg0_2.Draw()
	c0.SaveAs("results_phase%.4f_%d/invmass_D0Dbar0_fit_factor%s%.1f_%d.png"%(phase,event_n,type,event,event_n))

	c00 = ROOT.TCanvas("c00", "Bbar0 decay invmass", 1420,710); c00.Divide(2,1) 
	c00.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h00_1 = pltini("h00_1",bin,mKpi,  '\ m(K^+\pi^-)\\;[GeV/c^2]',     ROOT.kBlack)
	f00_1 = fit(min(mKpi),max(mKpi),4,0.89555,0.0473); h0_2.Fit(f0_2,"R+"); h00_1.Fit(f00_1,"R+")
	leg00_1 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg00_1.AddEntry(h00_1,"\ B^0 \\; data \\;(MC)", "lf"); leg00_1.AddEntry(f00_1,"\ BW \\; fit", "l"); legset(leg00_1)
	h00_1.SetMaximum(1.8*h00_1.GetMaximum())
	h00_1.Draw(); leg00_1.Draw() 

	c00.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h00_2 = pltini("h00_2",bin,mKpi_conj,  '\ m(K^-\pi^+)\\;[GeV/c^2]',     ROOT.kBlack)
	f00_2 = fit(min(mKpi_conj),max(mKpi_conj),4,0.89555,0.0473); h00_2.Fit(f00_2,"R+")
	leg00_2 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg00_2.AddEntry(h00_2,"\ \\bar{B}^0 \\; data \\;(MC)", "lf"); leg00_2.AddEntry(f00_2,"\ BW \\; fit", "l"); legset(leg00_2)
	h00_2.SetMaximum(h00_1.GetMaximum())
	h00_2.Draw(); leg00_2.Draw()
	c00.SaveAs("results_phase%.4f_%d/invmass_Kpi_fit_factor%s%.1f_%d.png"%(phase,event_n,type,event,event_n))

	with open('results_phase%.4f_%d/invmass_fit_parameters_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n),mode='a') as f:
		f.write("bin=%d A1=%d A2=%d A3=%d"%(bin,A1,A2,A3)+ "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(3770)",f0_1.GetParameters()[0],f0_1.GetParError(0), f0_1.GetParameters()[1],f0_1.GetParError(1),f0_1.GetChisquare(),"B0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(3770)",f0_2.GetParameters()[0],f0_2.GetParError(0), f0_2.GetParameters()[1],f0_2.GetParError(1),f0_2.GetChisquare(),"Bbar0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4040)",f0_1.GetParameters()[3],f0_1.GetParError(3), f0_1.GetParameters()[4],f0_1.GetParError(4),f0_1.GetChisquare(),"B0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4040)",f0_2.GetParameters()[3],f0_2.GetParError(3), f0_2.GetParameters()[4],f0_2.GetParError(4),f0_2.GetChisquare(),"Bbar0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4160)",f0_1.GetParameters()[6],f0_1.GetParError(6), f0_1.GetParameters()[7],f0_1.GetParError(7),f0_1.GetChisquare(),"B0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("psi(4160)",f0_2.GetParameters()[6],f0_2.GetParError(6), f0_2.GetParameters()[7],f0_2.GetParError(7),f0_2.GetChisquare(),"Bbar0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("K*0(892)", f00_1.GetParameters()[0],f00_1.GetParError(0), f00_1.GetParameters()[1],f00_1.GetParError(1),f00_1.GetChisquare(),"B0") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %.4f %s" % ("K*0(892)", f00_2.GetParameters()[0],f00_2.GetParError(0), f00_2.GetParameters()[1],f00_2.GetParError(1),f00_2.GetChisquare(),"Bbar0") + "\n")

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE        = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE       = float(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,PHASE,EVENT_N)