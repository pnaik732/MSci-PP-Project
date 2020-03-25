import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TH2F, TMath

name = "Dbarpi"; number = "13"

myFile1 = TFile.Open("Event_output.root", "READ")
hist1   = myFile1.Get("Data_s%s"%(number))

myFile2 = TFile.Open("Event_output.root", "READ")
hist2   = myFile2.Get("Model_cat1_s%s"%(number))

c1 = TCanvas("c1", "fit", 550,500)  
hist1.SetLineColor(ROOT.kWhite)
hist1.SetFillColorAlpha(ROOT.kOrange + 10, 0.4)
hist1.SetStats(False)

hist2.SetLineColor(ROOT.kViolet + 10)
hist2.SetLineWidth(2)
hist2.SetLineStyle(1)

hist2.SetMarkerStyle(20);
hist2.SetMarkerSize(0.5)
hist2.SetLineWidth(1);
hist2.SetMarkerColorAlpha(ROOT.kViolet + 10, 0.85);
hist2.SetStats(False)

#Create pad1
pad1 = TPad("pad1", "", 0, 0, 1, 1)
pad1.SetLeftMargin(0.1)
pad1.SetRightMargin(0.1)
pad1.SetTopMargin(0.1)
pad1.SetBottomMargin(0.1)
pad1.Draw()
pad1.cd()

hist1.Draw("B")
hist2.Draw("L same")

#Add the legend
legend = TLegend(0.55, 0.65, 0.85, 0.85)
legend.AddEntry(hist1, "LHCb data", "f")
legend.AddEntry(hist2, "amplitude fit", "lpe")

legend.SetTextFont(133)
legend.SetTextSize(20)
legend.SetLineColor(0)
legend.Draw()

c1.SaveAs("compare%s.png"%(name))

input()
