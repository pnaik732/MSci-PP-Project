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

def Fitting(variable, fit_func, weight, p0_list,bin):
    '''Fit the distribution with its analytic function, gives the fitting constant, used in Task1 and Task2.'''
    
    #---- Fit the histo distribution ----
    height, bin_edges = np.histogram (variable, bins=bin, weights=weight) #number distributions
    centre = bin_edges[:-1] + np.diff (bin_edges) / 2   #centre of each bin
    if len(p0_list) == 3:
        popt, pcov = curve_fit (fit_func, centre, height,p0=p0_list,bounds = ([0,0,0],[max(height),max(centre),1])) #finding constants and variance in the fitting function
    elif len(p0_list) == 9:
        popt, pcov = curve_fit (fit_func, centre, height,p0=p0_list,bounds = ([min(height),min(centre),0.00001,min(height),min(centre),0.00001,min(height),min(centre),0.00001],[max(height),max(centre),10,max(height),max(centre),0.1,max(height),max(centre),0.2])) #finding constants and variance in the fitting function

    #-------- plottings ------------------
    variable_fit = np.linspace (bin_edges[0], bin_edges[-1], 10000)  #array for independent variables
    return height, popt, np.sqrt(np.diag(pcov)), variable_fit

def find_yield(variable, variable_popt, variable_err, bin):
	delta_variable = (max(variable) - min(variable)) / bin
	yield_variable = (np.sqrt(2 * np.pi) * variable_popt[0] * variable_popt[2]) / delta_variable
	err1_yield     = np.sqrt(yield_variable) 
	err2_yield     = yield_variable * np.sqrt((variable_err[0] / variable_popt[0]) ** 2 + (variable_err[2] / variable_popt[2]) ** 2)
	return yield_variable, err1_yield, err2_yield
	    
def find_a_t(c_p, c_p_err1, c_p_err2, c_n, c_n_err1, c_n_err2):
	c_p_err = c_p_err1 + c_p_err2
	c_n_err = c_n_err1 + c_n_err2
	Asy     = (c_p-c_n) / (c_p+c_n)
	Asy_err = ((c_p_err * 2 * c_n / (c_p + c_n) ** 2) ** 2 + (c_n_err * 2 * c_p / (c_p + c_n) ** 2) ** 2) ** 0.5
	return Asy,Asy_err
	
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


	
def main (program,name):
	CM_data      = np.transpose(np.loadtxt('results_%s/B0_CM_variables.txt'%(name)))
	CM_data_conj = np.transpose(np.loadtxt('results_%s/B0_CM_variables_conj.txt'%(name)))

	mDDbar      = CM_data[1];            mDDbar_conj      = CM_data_conj[1]
	mKpi        = CM_data[2];            mKpi_conj        = CM_data_conj[2]
	C_T         = CM_data[5];            C_T_conj         = CM_data_conj[5]
	mDDbar_p    = CM_data[1][(C_T > 0)]; mDDbar_p_conj    = CM_data_conj[1][(C_T_conj > 0)]
	mKpi_p      = CM_data[2][(C_T > 0)]; mKpi_p_conj      = CM_data_conj[2][(C_T_conj > 0)]
	sWeight_p   = CM_data[6][(C_T > 0)]; sWeight_p_conj   = CM_data_conj[6][(C_T_conj > 0)] 
	mB0_p       = CM_data[7][(C_T > 0)]; mB0_p_conj       = CM_data_conj[7][(C_T_conj > 0)]
	mDDbar_n    = CM_data[1][(C_T < 0)]; mDDbar_n_conj    = CM_data_conj[1][(C_T_conj < 0)]
	mKpi_n      = CM_data[2][(C_T < 0)]; mKpi_n_conj      = CM_data_conj[2][(C_T_conj < 0)]
	sWeight_n   = CM_data[6][(C_T < 0)]; sWeight_n_conj   = CM_data_conj[6][(C_T_conj < 0)]
	mB0_n       = CM_data[7][(C_T < 0)]; mB0_n_conj       = CM_data_conj[7][(C_T_conj < 0)]
	
	bin=100
	f0,axs0 = plt.subplots(2,2,figsize=(9,8))
	p0_list = [1,1,1]
	(ax0_1, ax0_2), (ax0_3, ax0_4) = axs0
	if name ==  "run1ADDrun2_new_BW":
		function = BreitWigner ; fit_label = "fit (BW)"
	else:
		function = Gaussian;     fit_label = "fit (Gaussian)"
	
	mB0_p_height,      mB0_p_popt,      mB0_p_err      ,_ =  Fitting(np.array(mB0_p),      function, sWeight_p,      p0_list,bin)
	mB0_n_height,      mB0_n_popt,      mB0_n_err      ,_ =  Fitting(np.array(mB0_n),      function, sWeight_n,      p0_list,bin)
	mB0_p_height_conj, mB0_p_popt_conj, mB0_p_err_conj ,_ =  Fitting(np.array(mB0_p_conj), function, sWeight_p_conj, p0_list,bin)
	mB0_n_height_conj, mB0_n_popt_conj, mB0_n_err_conj ,_ =  Fitting(np.array(mB0_n_conj), function, sWeight_n_conj, p0_list,bin)

	mB0_p_yield,       mB0_p_err1_yield,      mB0_p_err2_yield      = find_yield(mB0_p,      mB0_p_popt,      mB0_p_err,      bin)
	mB0_n_yield,       mB0_n_err1_yield,      mB0_n_err2_yield      = find_yield(mB0_n,      mB0_n_popt,      mB0_n_err,      bin)
	mB0_p_yield_conj,  mB0_p_err1_yield_conj, mB0_p_err2_yield_conj = find_yield(mB0_p_conj, mB0_p_popt_conj, mB0_p_err_conj, bin)
	mB0_n_yield_conj,  mB0_n_err1_yield_conj, mB0_n_err2_yield_conj = find_yield(mB0_n_conj, mB0_n_popt_conj, mB0_n_err_conj, bin)

	A_T,      A_T_err     = find_a_t(mB0_p_yield,       mB0_p_err1_yield,      mB0_p_err2_yield,      mB0_n_yield,       mB0_n_err1_yield,      mB0_n_err2_yield)
	A_T_conj, A_T_err_conj = find_a_t(mB0_p_yield_conj,  mB0_p_err1_yield_conj, mB0_p_err2_yield_conj, mB0_n_yield_conj,  mB0_n_err1_yield_conj, mB0_n_err2_yield_conj)
	
	a_cp     = 0.5 * (A_T - A_T_conj)
	a_cp_err = 0.5 * np.sqrt(A_T_err ** 2 + A_T_err_conj ** 2)
	
	A_T_weights,      A_T_weights_err      = find_a_t(np.sum(sWeight_p),        np.sqrt(np.sum(sWeight_p)),           0,     np.sum(sWeight_n),        np.sqrt(np.sum(sWeight_n)),           0)
	A_T_weights_conj, A_T_weights_err_conj = find_a_t(np.sum(sWeight_p_conj),   np.sqrt(np.sum(sWeight_p_conj)),      0,     np.sum(sWeight_n_conj),   np.sqrt(np.sum(sWeight_n_conj)),      0)

	a_cp_weights     = 0.5 * (A_T_weights - A_T_weights_conj)
	a_cp_weights_err = 0.5 * np.sqrt(A_T_weights_err ** 2 + A_T_weights_err_conj ** 2)
	
	OneD_histo(ax0_1,np.array(mB0_p),      sWeight_p,      function, "$B^0$ ($C_T>0$)",                fit_label, "$m(D^0\\bar{D}^0K^+\pi^-)[GeV/c^2]$", "b-", "step", p0_list,bin)
	OneD_histo(ax0_2,np.array(mB0_n),      sWeight_n,      function, "$B^0$ ($C_T<0$)",                fit_label, "$m(D^0\\bar{D}^0K^+\pi^-)[GeV/c^2]$", "b-", "step", p0_list,bin)
	OneD_histo(ax0_3,np.array(mB0_p_conj), sWeight_p_conj, function, "$\\bar{B}^0$ ($-\\bar{C}_T>0$)", fit_label, "$m(\\bar{D}^0D^0K^-\pi^+)[GeV/c^2]$", "b-", "step", p0_list,bin)
	OneD_histo(ax0_4,np.array(mB0_n_conj), sWeight_n_conj, function, "$\\bar{B}^0$ ($-\\bar{C}_T<0$)", fit_label, "$m(\\bar{D}^0D^0K^-\pi^+)[GeV/c^2]$", "b-", "step", p0_list,bin)
	
	mB0_all_min        = [min(mB0_p),         min(mB0_n),        min(mB0_p_conj),        min(mB0_n_conj)]
	mB0_all_max        = [max(mB0_p),         max(mB0_n),        max(mB0_p_conj),        max(mB0_n_conj)]
	mB0_all_height_max = [max(mB0_p_height),  max(mB0_n_height), max(mB0_p_height_conj), max(mB0_n_height_conj)]

	plt.setp(axs0, xlim=(min(mB0_all_min),max(mB0_all_max)), ylim=(0,1.3*max(mB0_all_height_max)))
	f0.tight_layout()
	f0.savefig("results_%s/invmass_B0_fit_tp.png"%(name))
	
	event_n_weight = np.sum(mB0_p_height)+np.sum(mB0_n_height)+np.sum(mB0_p_height_conj)+np.sum(mB0_n_height_conj)
	print("total number of data after weight:", event_n_weight,"=",np.sum(mB0_p_height),"+",np.sum(mB0_n_height),"+",np.sum(mB0_p_height_conj),"+",np.sum(mB0_n_height_conj))
	print("weights", np.sum(sWeight_p),np.sum(sWeight_n),np.sum(sWeight_p_conj),np.sum(sWeight_n_conj))
	with open('results_%s/invmass_fit_tp_parameters.txt'%(name),mode='a') as f:
		f.write("bin=%d"%(bin)+ "\n")
		f.write("    A      err_A    mean   err_mean sigma  err_sigma type" +"\n")
		f.write("%s %.4f  %.4f  %.4f  %.4f   %.4f  %.4f    %s" % ("B0",       mB0_p_popt[0],        mB0_p_err[0],      mB0_p_popt[1],        mB0_p_err[1],         mB0_p_popt[2],        mB0_p_err[2],        "B0(C_T>0)")    + "\n")
		f.write("%s %.4f  %.4f  %.4f  %.4f   %.4f  %.4f    %s" % ("B0",       mB0_n_popt[0],        mB0_n_err[0],      mB0_n_popt[1],        mB0_n_err[1],         mB0_n_popt[2],        mB0_n_err[2],        "B0(C_T<0)")    + "\n")
		f.write("%s %.4f  %.4f  %.4f  %.4f   %.4f  %.4f    %s" % ("B0",       mB0_p_popt_conj[0],   mB0_p_err_conj[0], mB0_p_popt_conj[1],   mB0_p_err_conj[1],    mB0_p_popt_conj[2],   mB0_p_err_conj[2],   "Bbar0(-C_T>0)") + "\n")
		f.write("%s %.4f  %.4f  %.4f  %.4f   %.4f  %.4f    %s" % ("B0",       mB0_n_popt_conj[0],   mB0_n_err_conj[0], mB0_n_popt_conj[1],   mB0_n_err_conj[1],    mB0_n_popt_conj[2],   mB0_n_err_conj[2],   "Bbar0(-C_T<0)") + "\n")
		f.write("\n"+"=========================== Method 1 (from fittings) =================================="+"\n")
		f.write("\n"+"yields"+ "\n")
		f.write("events             eta           err_self     err_fit" + "\n")
		f.write("%s          %.4f          %.4f      %.4f"%("(C_T>0)",     mB0_p_yield,       mB0_p_err1_yield,      mB0_p_err2_yield)+ "\n")
		f.write("%s          %.4f          %.4f      %.4f"%("(C_T<0)",     mB0_n_yield,       mB0_n_err1_yield,      mB0_n_err2_yield)+ "\n")
		f.write("%s      %.4f          %.4f      %.4f"%("(-Cbar_T>0)", mB0_p_yield_conj,  mB0_p_err1_yield_conj, mB0_p_err2_yield_conj)+ "\n")
		f.write("%s      %.4f          %.4f      %.4f"%("(-Cbar_T<0)", mB0_n_yield_conj,  mB0_n_err1_yield_conj, mB0_n_err2_yield_conj)+ "\n")
		f.write("\n"+"Asymmetries"+"\n")
		f.write("A_T       errA_T    A_T_conj errA_T_conj   a_cp     err_a_cp"+ "\n")
		f.write("%.4f   %.4f   %.4f   %.4f       %.4f   %.4f"%(A_T,      A_T_err,A_T_conj, A_T_err_conj,a_cp,a_cp_err)+ "\n")
		f.write("\n"+"=========================== Method 2 (from weights) =================================="+"\n")
		f.write("\n" + "events      eta(weights)    err_self(weights) " + "\n")
		f.write("%s      %.4f         %.4f"%("(C_T>0)",     np.sum(sWeight_p),    np.sqrt(np.sum(sWeight_p)))+ "\n")
		f.write("%s      %.4f         %.4f"%("(C_T<0)",     np.sum(sWeight_n),    np.sqrt(np.sum(sWeight_n)))+ "\n")
		f.write("%s  %.4f         %.4f"%("(-Cbar_T>0)", np.sum(sWeight_p_conj),    np.sqrt(np.sum(sWeight_p_conj)))+ "\n")
		f.write("%s  %.4f         %.4f"%("(-Cbar_T<0)", np.sum(sWeight_n_conj),    np.sqrt(np.sum(sWeight_n_conj)))+ "\n")
		f.write("\n"+"Asymmetries (using weights)"+"\n")
		f.write("A_T       errA_T    A_T_conj errA_T_conj   a_cp     err_a_cp"+ "\n")
		f.write("%.4f   %.4f   %.4f   %.4f       %.4f   %.4f"%(A_T_weights,      A_T_weights_err, A_T_weights_conj, A_T_weights_err_conj,a_cp_weights,a_cp_weights_err)+ "\n")

	
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    NAME        = str(sys.argv[1])
    main(PROGNAME,NAME)
