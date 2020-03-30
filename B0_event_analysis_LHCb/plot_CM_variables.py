#************************************************
#plot the distributions of CM variables

#@Author: Tailin Zhu
#************************************************

import numpy as np
import ROOT
from ROOT import gPad
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys

### generate plots and add data ###
def pltini(no,bin,data,sweight,xaxis,color):
    h = ROOT.TH1F(no, "", bin, min(data), max(data))
    for m in range (len(data)):
	    h.Fill(data[m],sweight[m])
    pltset (h,xaxis,color)
    return h

### customise the front size and colour ###    
def pltset (input, xaxis,color):
	#function to set up the fontsize of the axis labels and titles.
	gr = input
	labelsize = 0.05
	gr.SetStats(False)
	gr.SetLineColor(color)
	if color == ROOT.kRed:
		gr.SetLineStyle(2)
	gr.GetXaxis().SetTitle(xaxis)
	gr.GetYaxis().SetTitle('Candidates')
	gr.GetXaxis().SetLabelFont(132)
	gr.GetYaxis().SetLabelFont(132)
	gr.GetXaxis().SetTitleFont(132)
	gr.GetYaxis().SetTitleFont(132)
	gr.GetXaxis().SetLabelSize(labelsize)
	gr.GetXaxis().SetTitleSize(labelsize)
	gr.GetYaxis().SetLabelSize(labelsize)
	gr.GetYaxis().SetTitleSize(labelsize)
	gr.SetMinimum(0)

### add legend to the plot ###
def legset(input):
	legend = input
	size = 25
	legend.SetTextFont(133)
	legend.SetTextSize(size)
	legend.SetLineColor(0)
    
def main(program,name):
	#0.phi, 1. mddbar, 2. mkpi, 3. costhetaddbar, 4. costhetakpi, 5.C_T
	CM_data      = np.transpose(np.loadtxt('results_%s/B0_CM_variables.txt'%(name)))
	CM_data_conj = np.transpose(np.loadtxt('results_%s/B0_CM_variables_conj.txt'%(name)))

	mDDbar      = CM_data[1];            mDDbar_conj      = CM_data_conj[1]
	mKpi        = CM_data[2];            mKpi_conj        = CM_data_conj[2]
	C_T         = CM_data[5];            C_T_conj         = CM_data_conj[5]
	phi_p       = CM_data[0][(C_T > 0)]; phi_p_conj       = CM_data_conj[0][(C_T_conj > 0)]
	mDDbar_p    = CM_data[1][(C_T > 0)]; mDDbar_p_conj    = CM_data_conj[1][(C_T_conj > 0)]
	mKpi_p      = CM_data[2][(C_T > 0)]; mKpi_p_conj      = CM_data_conj[2][(C_T_conj > 0)]
	costhetaD_p = CM_data[3][(C_T > 0)]; costhetaD_p_conj = CM_data_conj[3][(C_T_conj > 0)]
	costhetaK_p = CM_data[4][(C_T > 0)]; costhetaK_p_conj = CM_data_conj[4][(C_T_conj > 0)]
	sWeight_p   = CM_data[6][(C_T > 0)]; sWeight_p_conj   = CM_data_conj[6][(C_T_conj > 0)] 
	mB0_p       = CM_data[7][(C_T > 0)]; mB0_p_conj       = CM_data_conj[7][(C_T_conj > 0)]
	phi_n       = CM_data[0][(C_T < 0)]; phi_n_conj       = CM_data_conj[0][(C_T_conj < 0)]
	mDDbar_n    = CM_data[1][(C_T < 0)]; mDDbar_n_conj    = CM_data_conj[1][(C_T_conj < 0)]
	mKpi_n      = CM_data[2][(C_T < 0)]; mKpi_n_conj      = CM_data_conj[2][(C_T_conj < 0)]
	costhetaD_n = CM_data[3][(C_T < 0)]; costhetaD_n_conj = CM_data_conj[3][(C_T_conj < 0)]
	costhetaK_n = CM_data[4][(C_T < 0)]; costhetaK_n_conj = CM_data_conj[4][(C_T_conj < 0)]
	sWeight_n   = CM_data[6][(C_T < 0)]; sWeight_n_conj   = CM_data_conj[6][(C_T_conj < 0)]
	mB0_n       = CM_data[7][(C_T < 0)]; mB0_n_conj       = CM_data_conj[7][(C_T_conj < 0)]
	
	bin = 100

	#*****************************************************************************************
	#invariant mass  -  m(D0 Dbar0)
	#*****************************************************************************************
	c1 = ROOT.TCanvas("c1", "B0 decay invmass", 1420,710); c1.Divide(2,1) 
	c1.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h1_1 = pltini("h1_1",bin,mDDbar_p,sWeight_p,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kAzure+2)
	h1_2 = pltini("h1_2",bin,mDDbar_n,sWeight_n,'\ m(D^0\\bar{D}^0)\\;[GeV/c^2]',ROOT.kRed)
	leg1 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg1.AddEntry(h1_1,"\ B^0(C_T > 0)", "lf"); leg1.AddEntry(h1_2,"\ B^0(C_T < 0)", "lf"); legset(leg1)
	h1_1.SetMaximum(1.8*h1_1.GetMaximum())
	h1_1.Draw("HIST"); h1_2.Draw("HIST same"); leg1.Draw()

	c1.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h11_1 = pltini("h11_1",bin,mDDbar_p_conj,sWeight_p_conj,'\ m(\\bar{D}^0D^0)\\;[GeV/c^2]',ROOT.kAzure+2)
	h11_2 = pltini("h11_2",bin,mDDbar_n_conj,sWeight_n_conj,'\ m(\\bar{D}^0D^0)\\;[GeV/c^2]',ROOT.kRed)
	leg11 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg11.AddEntry(h11_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg11.AddEntry(h11_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg11)
	h11_1.SetMaximum(h1_1.GetMaximum())
	h11_1.Draw("HIST"); h11_2.Draw("HIST same"); leg11.Draw()
	c1.SaveAs("results_%s/invmass_D0Dbar0.png"%(name))

	#*****************************************************************************************
	#invariant mass  -  m(K pi)
	#*****************************************************************************************
	c2 = ROOT.TCanvas("c2", "B0 decay invmass", 1420,710); c2.Divide(2,1) 
	c2.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h2_1 = pltini("h2_1",bin,mKpi_p,sWeight_p,'\ m(K^+\pi^-)\\;[GeV/c^2]',ROOT.kAzure+2)
	h2_2 = pltini("h2_2",bin,mKpi_n,sWeight_n,'\ m(K^+\pi^-)\\;[GeV/c^2]',ROOT.kRed)
	leg2 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg2.AddEntry(h2_1,"\ B^0(C_T > 0)", "lf"); leg2.AddEntry(h2_2,"\ B^0(C_T < 0)", "lf"); legset(leg2)
	h2_1.SetMaximum(1.8*h2_1.GetMaximum())
	h2_1.Draw("HIST"); h2_2.Draw("HIST same"); leg2.Draw()

	c2.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h22_1 = pltini("h22_1",bin,mKpi_p_conj,sWeight_p_conj,'\ m(K^-\pi^+)\\;[GeV/c^2]',ROOT.kAzure+2)
	h22_2 = pltini("h22_2",bin,mKpi_n_conj,sWeight_n_conj,'\ m(K^-\pi^+)\\;[GeV/c^2]',ROOT.kRed)
	leg22 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg22.AddEntry(h22_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg22.AddEntry(h22_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg22)
	h22_1.SetMaximum(h2_1.GetMaximum())
	h22_1.Draw("HIST"); h22_2.Draw("HIST same"); leg22.Draw()
	c2.SaveAs("results_%s/invmass_KpPim.png"%(name))

	#*****************************************************************************************
	#helicity angle  -  cos(theta_D)
	#*****************************************************************************************
	c3 = ROOT.TCanvas("c3", "B0 decay helicity anlge", 1420,710); c3.Divide(2,1) 
	c3.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h3_1 = pltini("h3_1",bin,costhetaD_p,sWeight_p,'\ \cos(\\theta_D)',ROOT.kAzure+2)
	h3_2 = pltini("h3_2",bin,costhetaD_n,sWeight_n,'\ \cos(\\theta_D)',ROOT.kRed)
	leg3 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg3.AddEntry(h3_1,"\ B^0(C_T > 0)", "lf"); leg3.AddEntry(h3_2,"\ B^0(C_T < 0)", "lf"); legset(leg3)
	h3_1.SetMaximum(1.8*h3_1.GetMaximum())
	h3_1.Draw("HIST"); h3_2.Draw("HIST same"); leg3.Draw()

	c3.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h33_1 = pltini("h33_1",bin,costhetaD_p_conj,sWeight_p_conj,'\ \cos(\\theta_D)',ROOT.kAzure+2)
	h33_2 = pltini("h33_2",bin,costhetaD_n_conj,sWeight_n_conj,'\ \cos(\\theta_D)',ROOT.kRed)
	leg33 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg33.AddEntry(h33_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg33.AddEntry(h33_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg33)
	h33_1.SetMaximum(h3_1.GetMaximum())
	h33_1.Draw("HIST"); h33_2.Draw("HIST same"); leg33.Draw()
	c3.SaveAs("results_%s/helangle_D0Dbar0.png"%(name))

	#*****************************************************************************************
	#helicity angle  -  cos(theta_K)
	#*****************************************************************************************
	c4 = ROOT.TCanvas("c4", "B0 decay helicity anlge", 1420,710); c4.Divide(2,1) 
	c4.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h4_1 = pltini("h4_1",bin,costhetaK_p,sWeight_p,'\ \cos(\\theta_K)',ROOT.kAzure+2)
	h4_2 = pltini("h4_2",bin,costhetaK_n,sWeight_n,'\ \cos(\\theta_K)',ROOT.kRed)
	leg4 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg4.AddEntry(h4_1,"\ B^0(C_T > 0)", "lf"); leg4.AddEntry(h4_2,"\ B^0(C_T < 0)", "lf"); legset(leg4)
	h4_1.SetMaximum(1.8*h4_1.GetMaximum())
	h4_1.Draw("HIST"); h4_2.Draw("HIST same"); leg4.Draw()

	c4.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h44_1 = pltini("h44_1",bin,costhetaK_p_conj,sWeight_p_conj,'\ \cos(\\theta_K)',ROOT.kAzure+2)
	h44_2 = pltini("h44_2",bin,costhetaK_n_conj,sWeight_n_conj,'\ \cos(\\theta_K)',ROOT.kRed)
	leg44 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg44.AddEntry(h44_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg44.AddEntry(h44_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg44)
	h44_1.SetMaximum(h4_1.GetMaximum())
	h44_1.Draw("HIST"); h44_2.Draw("HIST same"); leg44.Draw()
	c4.SaveAs("results_%s/helangle_KpPim.png"%(name))

	#*****************************************************************************************
	#plane angle  -  phi
	#*****************************************************************************************
	c5 = ROOT.TCanvas("c5", "B0 decay plane angle", 1420,710); c5.Divide(2,1) 
	c5.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h5_1 = pltini("h5_1",bin,phi_p,sWeight_p,'\ \phi\\;[rad]',ROOT.kAzure+2)
	h5_2 = pltini("h5_2",bin,phi_n,sWeight_n,'\ \phi\\;[rad]',ROOT.kRed)
	leg5 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg5.AddEntry(h5_1,"\ B^0(C_T > 0)", "lf"); leg5.AddEntry(h5_2,"\ B^0(C_T < 0)", "lf"); legset(leg5)
	h5_1.SetMaximum(1.8*h5_1.GetMaximum())
	h5_1.Draw("HIST"); h5_2.Draw("HIST same"); leg5.Draw()

	c5.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h55_1 = pltini("h55_1",bin,phi_p_conj,sWeight_p_conj,'\ \phi\\;[rad]',ROOT.kAzure+2)
	h55_2 = pltini("h55_2",bin,phi_n_conj,sWeight_n_conj,'\ \phi\\;[rad]',ROOT.kRed)
	leg55 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg55.AddEntry(h55_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg55.AddEntry(h55_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg55)
	h55_1.SetMaximum(h5_1.GetMaximum())
	h55_1.Draw("HIST"); h55_2.Draw("HIST same"); leg55.Draw()
	c5.SaveAs("results_%s/planeangle.png"%(name))
	
	#*****************************************************************************************
	#invariant mass  -  m(D0 Dbar0 K pi)
	#*****************************************************************************************
	c6 = ROOT.TCanvas("c6", "B0 decay invmass", 1420,710); c6.Divide(2,1) 
	c6.cd(1); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h6_1 = pltini("h6_1",bin,mB0_p,sWeight_p,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kAzure+2)
	h6_2 = pltini("h6_2",bin,mB0_n,sWeight_n,'\ m(D^0\\bar{D}^0K^+\pi^-)\\;[GeV/c^2]',ROOT.kRed)
	leg6 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg6.AddEntry(h6_1,"\ B^0(C_T > 0)", "lf"); leg6.AddEntry(h6_2,"\ B^0(C_T < 0)", "lf"); legset(leg6)
	h6_1.SetMaximum(1.8*h6_1.GetMaximum())
	h6_1.Draw("HIST"); h6_2.Draw("HIST same"); leg6.Draw()

	c6.cd(2); gPad.SetMargin(0.2,0.05,0.15,0.1)
	h66_1 = pltini("h66_1",bin,mB0_p_conj,sWeight_p_conj,'\ m(\\bar{D}^0D^0K^-\pi^+)\\;[GeV/c^2]',ROOT.kAzure+2)
	h66_2 = pltini("h66_2",bin,mB0_n_conj,sWeight_n_conj,'\ m(\\bar{D}^0D^0K^-\pi^+)\\;[GeV/c^2]',ROOT.kRed)
	leg66 = ROOT.TLegend(0.6, 0.74, 0.89, 0.89); leg66.AddEntry(h66_1,"\ \\bar{B}^0(-\\bar{C}_T > 0)", "lf"); leg66.AddEntry(h66_2,"\ \\bar{B}^0(-\\bar{C}_T < 0)", "lf"); legset(leg66)
	h66_1.SetMaximum(h6_1.GetMaximum())
	h66_1.Draw("HIST"); h66_2.Draw("HIST same"); leg66.Draw()
	c6.SaveAs("results_%s/invmass_B0.png"%(name))


	#combined figure
	img1 = mpimg.imread("results_%s/invmass_D0Dbar0.png"%(name))
	img2 = mpimg.imread("results_%s/invmass_KpPim.png"%(name))
	img3 = mpimg.imread("results_%s/helangle_D0Dbar0.png"%(name))
	img4 = mpimg.imread("results_%s/helangle_KpPim.png"%(name))
	img5 = mpimg.imread("results_%s/planeangle.png"%(name))

	fig1,axs1 = plt.subplots(2,1,figsize=(60,50))
	(ax1),(ax2) = axs1

	fig2,axs2 = plt.subplots(3,1,figsize=(60,75))
	(ax3),(ax4),(ax5) = axs2

	fig3,axs3 = plt.subplots(5,1,figsize=(60,125))
	(ax1_1),(ax2_1),(ax3_1),(ax4_1),(ax5_1) = axs3

	ax1.imshow(img1); ax1.set_axis_off()
	ax2.imshow(img2); ax2.set_axis_off()
	ax3.imshow(img3); ax3.set_axis_off()
	ax4.imshow(img4); ax4.set_axis_off()
	ax5.imshow(img5); ax5.set_axis_off()

	ax1_1.imshow(img1); ax1_1.set_axis_off()
	ax2_1.imshow(img2); ax2_1.set_axis_off()
	ax3_1.imshow(img3); ax3_1.set_axis_off()
	ax4_1.imshow(img4); ax4_1.set_axis_off()
	ax5_1.imshow(img5); ax5_1.set_axis_off()

	fig1.tight_layout()
	fig2.tight_layout()
	fig3.tight_layout()
	fig1.savefig("results_%s/CM_variable_invariant_mass.png"%(name))
	fig2.savefig("results_%s/CM_variable_angles.png"%(name))
	fig3.savefig("results_%s/CM_variable_all.png"%(name))

	#input()        #python3
	#raw_input()    #python2
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    NAME        = str(sys.argv[1])
    main(PROGNAME,NAME)
