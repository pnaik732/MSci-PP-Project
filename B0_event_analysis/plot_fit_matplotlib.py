import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import sys

def Gaussian(x, A, mean, sigma):
    '''The analytical function of Gaussian distribution.'''
    return A*np.exp(-(((x-mean)/sigma)**2)/2)

def BreitWigner (E, A, M, Gamma):
    gamma = (M**2*(M**2+Gamma**2))**0.5
    k = (2**1.8*M*Gamma*gamma)/(np.pi*(M**2+gamma)**0.5)
    f = A * k/((E**2-M**2)**2+(M*Gamma)**2)
    return f
    
def BreitWigner3 (E, A1, M1, Gamma1,A2, M2, Gamma2,A3, M3, Gamma3):
    f1 = BreitWigner (E, A1, M1, Gamma1)
    f2 = BreitWigner (E, A2, M2, Gamma2)
    f3 = BreitWigner (E, A3, M3, Gamma3)
    f = f1 + f2 + f3
    return f

def Gaussian3 (x, A1, mean1, sigma1,A2, mean2, sigma2,A3, mean3, sigma3):
    f1 = Gaussian (x, A1, mean1, sigma1)
    f2 = Gaussian (x, A2, mean2, sigma2)
    f3 = Gaussian (x, A3, mean3, sigma3)
    f = f1 + f2 + f3
    return f

def Fitting(variable, fit_func, weight, p0_list,bin):
    '''Fit the distribution with its analytic function, gives the fitting constant, used in Task1 and Task2.'''
    
    #---- Fit the histo distribution ----
    height, bin_edges = np.histogram (variable, bins=bin, weights=weight) #number distributions
    centre = bin_edges[:-1] + np.diff (bin_edges) / 2   #centre of each bin
    if len(p0_list) == 3:
        popt, pcov = curve_fit (fit_func, centre, height,p0=p0_list,bounds = ([min(height),min(centre),0.00001],[max(height),max(centre),0.2])) #finding constants and variance in the fitting function
    elif len(p0_list) == 9:
        popt, pcov = curve_fit (fit_func, centre, height,p0=p0_list,bounds = ([min(height),min(centre),0.00001,min(height),min(centre),0.00001,min(height),min(centre),0.00001],[max(height),max(centre),10,max(height),max(centre),0.1,max(height),max(centre),0.2])) #finding constants and variance in the fitting function

    #-------- plottings ------------------
    variable_fit = np.linspace (bin_edges[0], bin_edges[-1], 10000)  #array for independent variables
    return height, popt, np.sqrt(np.diag(pcov)), variable_fit
    
def OneD_histo(ax,variable, weight, fit_func, title, label0, xlabel, style, histtype, p0_list, bin):
	'''Plot the histogram for the interested variable, used in Task1 and Task2.'''
	_, popt, perr, variable_fit = Fitting(variable, fit_func,weight, p0_list,bin)
	nstd = 1. # to draw 5sigma intervals
	popt_up = popt + nstd * perr
	popt_dw = popt - nstd * perr

	fit = fit_func(variable_fit, *popt)
	fit_up = fit_func(variable_fit, *popt_up)
	fit_dw = fit_func(variable_fit, *popt_dw)

	ax.hist(variable, bins=bin, histtype=histtype, label=title,lw=1, weights=weight, color="k")  #histo distribution
	ax.plot(variable_fit, fit_func(variable_fit, *popt), style, label = label0,linewidth = 1) #fitting curve
	ax.fill_between(variable_fit, fit_up, fit_dw, alpha=0.25, color="b", label='$1\sigma$ interval')
	ax.set(xlabel=xlabel, ylabel="Candidates")
	ax.xaxis.set_label_coords(0.8, -0.11,  transform=ax.transAxes)
	ax.yaxis.set_label_coords(-0.12, 0.9, transform=ax.transAxes)
	ax.legend(fontsize=13,frameon=False)

def main (program,type,event,phase,event_n):
	CM_data      = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n)))
	CM_data_conj = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d_conj.txt'%(phase,event_n,type,event,event_n)))

	mDDbar      = CM_data[1];            mDDbar_conj      = CM_data_conj[1]
	mKpi        = CM_data[2];            mKpi_conj        = CM_data_conj[2]
	C_T         = CM_data[5];            C_T_conj         = CM_data_conj[5]
	mDDbar_p    = CM_data[1][(C_T > 0)]; mDDbar_p_conj    = CM_data_conj[1][(C_T_conj > 0)]
	mKpi_p      = CM_data[2][(C_T > 0)]; mKpi_p_conj      = CM_data_conj[2][(C_T_conj > 0)]
	mDDbar_n    = CM_data[1][(C_T < 0)]; mDDbar_n_conj    = CM_data_conj[1][(C_T_conj < 0)]
	mKpi_n      = CM_data[2][(C_T < 0)]; mKpi_n_conj      = CM_data_conj[2][(C_T_conj < 0)]

	bin=100

	f1,axs1 = plt.subplots(2,2,figsize=(9,8))
	(ax1_1, ax1_2), (ax1_3, ax1_4) = axs1
	p1_list = [50,0.89555,0.0473]
	mKpi_p_p_height,      mKpi_p_popt,      mKpi_p_err      ,_ =  Fitting(np.array(mKpi_p),      BreitWigner, np.ones(len(mKpi_p)),      p1_list,bin)
	mKpi_p_n_height,      mKpi_n_popt,      mKpi_n_err      ,_ =  Fitting(np.array(mKpi_n),      BreitWigner, np.ones(len(mKpi_n)),      p1_list,bin)
	mKpi_p_p_height_conj, mKpi_p_popt_conj, mKpi_p_err_conj ,_ =  Fitting(np.array(mKpi_p_conj), BreitWigner, np.ones(len(mKpi_p_conj)), p1_list,bin)
	mKpi_p_n_height_conj, mKpi_n_popt_conj, mKpi_n_err_conj ,_ =  Fitting(np.array(mKpi_n_conj), BreitWigner, np.ones(len(mKpi_n_conj)), p1_list,bin)
	OneD_histo(ax1_1,np.array(mKpi_p),      np.ones(len(mKpi_p)),      BreitWigner, "$B^0$ ($C_T>0$)",        "BW fit", "$m(K^+\pi^-)[GeV/c^2]$", "b-", "step",p1_list,bin)
	OneD_histo(ax1_2,np.array(mKpi_n),      np.ones(len(mKpi_n)),      BreitWigner, "$B^0$ ($C_T<0$)",        "BW fit", "$m(K^+\pi^-)[GeV/c^2]$", "b-", "step",p1_list,bin)
	OneD_histo(ax1_3,np.array(mKpi_p_conj), np.ones(len(mKpi_p_conj)), BreitWigner, "$\\bar{B}_0$ (-\\bar{C}_T>0$)", "BW fit", "$m(K^-\pi^+)[GeV/c^2]$", "b-", "step",p1_list,bin)
	OneD_histo(ax1_4,np.array(mKpi_n_conj), np.ones(len(mKpi_n_conj)), BreitWigner, "$\\bar{B}_0$ (-\\bar{C}_T<0$)", "BW fit", "$m(K^-\pi^+)[GeV/c^2]$", "b-", "step",p1_list,bin)
	plt.setp(axs1, xlim=(min(mKpi_p),max(mKpi_p)), ylim=(0,2*max(mKpi_p_height)))
	f1.tight_layout()
	f1.savefig("results_phase%.4f_%d/invmass_Kpi_fit_tp_factor%s%.1f_%d.png"%(phase,event_n,type,event,event_n))

	f2,axs2 = plt.subplots(2,2,figsize=(9,8))
	(ax2_1, ax2_2), (ax2_3, ax2_4) = axs2
	p2_list = [50,3.770,0.0272,10,4.039,0.080,10,4.191,0.070]
	mDDbar_p_height,      mDDbar_p_popt,      mDDbar_p_err      ,_ =  Fitting(np.array(mDDbar_p),      BreitWigner3, np.ones(len(mDDbar_p)),     p2_list,bin)
	mDDbar_n_height,      mDDbar_n_popt,      mDDbar_n_err      ,_ =  Fitting(np.array(mDDbar_n),      BreitWigner3, np.ones(len(mDDbar_n)),     p2_list,bin)
	mDDbar_p_height_conj, mDDbar_p_popt_conj, mDDbar_p_err_conj ,_ =  Fitting(np.array(mDDbar_p_conj), BreitWigner3, np.ones(len(mDDbar_p_conj)),p2_list,bin)
	mDDbar_n_height_conj, mDDbar_n_popt_conj, mDDbar_n_err_conj ,_ =  Fitting(np.array(mDDbar_n_conj), BreitWigner3, np.ones(len(mDDbar_n_conj)),p2_list,bin)
	OneD_histo(ax2_1,np.array(mDDbar_p),      np.ones(len(mDDbar_p)),      BreitWigner3, "$B^0$ ($C_T>0$)",        "BW fit", "$m(D^0\\bar{D}^0)[GeV/c^2]$", "b-", "step",p2_list,bin)
	OneD_histo(ax2_2,np.array(mDDbar_n),      np.ones(len(mDDbar_n)),      BreitWigner3, "$B^0$ ($C_T<0$)",        "BW fit", "$m(D^0\\bar{D}^0)[GeV/c^2]$", "b-", "step",p2_list,bin)
	OneD_histo(ax2_3,np.array(mDDbar_p_conj), np.ones(len(mDDbar_p_conj)), BreitWigner3, "$\\bar{B}_0$ (-\\bar{C}_T>0$)", "BW fit", "$m(\\bar{D}^0D^0)[GeV/c^2]$", "b-", "step",p2_list,bin)
	OneD_histo(ax2_4,np.array(mDDbar_n_conj), np.ones(len(mDDbar_n_conj)), BreitWigner3, "$\\bar{B}_0$ (-\\bar{C}_T<0$)", "BW fit", "$m(\\bar{D}^0D^0)[GeV/c^2]$", "b-", "step",p2_list,bin)
	plt.setp(axs2, xlim=(min(mDDbar_p),max(mDDbar_p)), ylim=(0,2*max(mDDbar_p)))
	f2.tight_layout()
	f2.savefig("results_phase%.4f_%d/invmass_D0Dbar0_fit_tp_factor%s%.1f_%d.png"%(phase,event_n,type,event,event_n))
    
	with open('results_phase%.4f_%d/invmass_fit_tp_parameters_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n),mode='a') as f:
		f.write("bin=%d"%(bin)+ "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(3770)",mDDbar_p_popt[1],     mDDbar_p_err[1],      mDDbar_p_popt[2],     mDDbar_p_err[2],     "B0(C_T>0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(3770)",mDDbar_n_popt[1],     mDDbar_n_err[1],      mDDbar_n_popt[2],     mDDbar_n_err[2],     "B0(C_T<0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(3770)",mDDbar_p_popt_conj[1],mDDbar_p_err_conj[1], mDDbar_p_popt_conj[2],mDDbar_p_err_conj[2],"Bbar0(C_T>0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(3770)",mDDbar_n_popt_conj[1],mDDbar_n_err_conj[1], mDDbar_n_popt_conj[2],mDDbar_n_err_conj[2],"Bbar0(C_T<0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4040)",mDDbar_p_popt[4],     mDDbar_p_err[4],      mDDbar_p_popt[5],     mDDbar_p_err[5],     "B0(C_T>0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4040)",mDDbar_n_popt[4],     mDDbar_n_err[4],      mDDbar_n_popt[5],     mDDbar_n_err[5],     "B0(C_T<0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4040)",mDDbar_p_popt_conj[4],mDDbar_p_err_conj[4], mDDbar_p_popt_conj[5],mDDbar_p_err_conj[5],"Bbar0(C_T>0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4040)",mDDbar_n_popt_conj[4],mDDbar_n_err_conj[4], mDDbar_n_popt_conj[5],mDDbar_n_err_conj[5],"Bbar0(C_T<0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4160)",mDDbar_p_popt[7],     mDDbar_p_err[7],      mDDbar_p_popt[8],     mDDbar_p_err[8],     "B0(C_T>0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4160)",mDDbar_n_popt[7],     mDDbar_n_err[7],      mDDbar_n_popt[8],     mDDbar_n_err[8],     "B0(C_T<0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4160)",mDDbar_p_popt_conj[7],mDDbar_p_err_conj[7], mDDbar_p_popt_conj[8],mDDbar_p_err_conj[8],"Bbar0(C_T>0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("psi(4160)",mDDbar_n_popt_conj[7],mDDbar_n_err_conj[7], mDDbar_n_popt_conj[8],mDDbar_n_err_conj[8],"Bbar0(C_T<0)") + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("K*0(892)", mKpi_p_popt[1],       mKpi_p_err[1],        mKpi_p_popt[2],       mKpi_p_err[2],       "B0(C_T>0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("K*0(892)", mKpi_n_popt[1],       mKpi_n_err[1],        mKpi_n_popt[2],       mKpi_n_err[2],       "B0(C_T<0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("K*0(892)", mKpi_p_popt_conj[1],  mKpi_p_err_conj[1],   mKpi_p_popt_conj[2],  mKpi_p_err_conj[2],  "B0(C_T>0)")    + "\n")
		f.write("%s %.4f %.4f %.4f %.4f %s" % ("K*0(892)", mKpi_n_popt_conj[1],  mKpi_n_err_conj[1],   mKpi_n_popt_conj[2],  mKpi_n_err_conj[2],  "B0(C_T>0)")    + "\n")

	# kde = stats.gaussian_kde(mKpi_p)
	# df = pd.DataFrame(mKpi_p)
	# ax1 = df.plot.kde()
	# xs = np.linspace(min(mKpi_p), max(mKpi_p), num=50)
	# y1 = kde.evaluate(xs)
	# fig, ax = plt.subplots()
	# ax.hist(mKpi_p,bins=100, histtype = "step", lw=1, label="$B_0$ ($C_T>0$)", weights = np.ones(len(mDDbar_p)),color = "k",normed=1)
	# ax.plot(xs, y1, label='Scott (default)')
	# ax.legend()
	
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE        = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE       = float(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,PHASE,EVENT_N)
