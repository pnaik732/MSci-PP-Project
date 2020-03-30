import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TH2F, TMath
import sys

def main(progname, file_name, number):
	myFile1 = TFile.Open("%s"%(file_name), "READ")
	hist1   = myFile1.Get("Data_s%s"%(number))

	myFile2 = TFile.Open("%s"%(file_name), "READ")
	hist2   = myFile2.Get("Model_cat1_s%s"%(number))

	c1 = TCanvas("c1", "fit", 1000,4000) 
	hist1.Sumw2()
	hist1.SetLineColor(ROOT.kOrange + 10)
	hist1.SetMarkerStyle(20)
	hist1.SetLineWidth(1)
	hist1.SetMarkerSize(0.7)
	hist1.SetMarkerColorAlpha(ROOT.kOrange + 10, 1);
	# hist1.SetFillColorAlpha(ROOT.kOrange + 10, 0.4)
	hist1.SetStats(False)


	hist2.SetLineColor(ROOT.kViolet + 10)
	hist2.SetLineStyle(1)
	hist2.Sumw2()
	hist2.SetMarkerStyle(20);
	hist2.SetMarkerSize(0.7)
	hist2.SetLineWidth(1);
	hist2.SetMarkerColorAlpha(ROOT.kViolet + 10, 1);
	hist2.SetStats(False)

	#Create pad1
	pad1 = TPad("pad1", "", 0, 0, 1, 1)
	pad1.SetLeftMargin(0.1)
	pad1.SetRightMargin(0.1)
	pad1.SetTopMargin(0.1)
	pad1.SetBottomMargin(0.1)
	pad1.Draw()
	pad1.cd()

	hist1.Draw("E1")
	hist2.Draw("E1 same")

	#Add the legend
	legend = TLegend(0.65, 0.75, 0.89, 0.89)
	legend.AddEntry(hist1, "LHCb data", "lpe")
	legend.AddEntry(hist2, "AmpGen fitter", "lpe")

	legend.SetTextFont(133)
	legend.SetTextSize(30)
	legend.SetLineColor(0)
	legend.Draw()

	c1.SaveAs("compare%s.png"%(number))
	
if __name__ == '__main__':
    PROGNAME    =     sys.argv[0]
    FILE_NAME   = str(sys.argv[1])
    NUMBER      = str(sys.argv[2])      
    main(PROGNAME,FILE_NAME,NUMBER)