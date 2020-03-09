#************************************************
#Binned Phase-Space Analysis for Four-body Decays 
#Author: Tailin Zhu
#Date: 6 March 2020

#Use 5 CM variables
#Split into 32 phase-space regions
#************************************************

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ROOT
import sys

### construct the bins containing equal number of events ###
def binning_patches(x,nbin):
    npt       = len(x)
    bin_edges = np.interp(np.linspace(0, npt, nbin + 1),
                          np.arange(npt),
                          np.sort(x))
    x_bins, x_patches = np.histogram(x, bin_edges)
    return x_patches

### binning the next variable from the events allocated in different bins of the previous variable ###   
def binning_data(x, y, x_patches,nbin,num):
    x_index  = np.digitize(x, x_patches[:-1])  #gives different number for event in different bins
    y_binned = []
    for i in range (1, nbin+1):
        index,    = np.where(x_index == i)    #find the index of where bin #1 or #2 occurs.
        y_partial = np.array(y[num+index])    
        y_binned.append(y_partial)
    y_binned = np.array(y_binned)
    return y_binned

### some repeated code used in find_bins() ###
def binning_data_more(y, x_binned, nbin, event_num,i):
	num = 0
	x_patches_list = []
	y_binned_total = []
	for bin in range(len(x_binned)):
		x_patches = binning_patches(x_binned[bin],nbin)
		y_binned  = binning_data(x_binned[bin],y,x_patches,nbin,num)
		x_patches_list.append(x_patches)
		y_binned_total.append(y_binned)
		num += int(event_num/(nbin**i))
	y_binned_total = np.concatenate(y_binned_total)
	return x_patches_list, y_binned_total

### binning each variable lists ###
def find_bins(data):
	nbin = 2                  #each step split into two bins
	event_num = len(data)     #total number of events
	data = np.transpose(data)

	x1 = data[0]
	x2 = data[1]
	x3 = data[2]
	x4 = data[3]
	x5 = data[4]

	x1_binned_total = np.array([x1])
	x1_patches_list, x2_binned_total = binning_data_more(x2, x1_binned_total, nbin, event_num,0)
	x2_patches_list, x3_binned_total = binning_data_more(x3, x2_binned_total, nbin, event_num,1)
	x3_patches_list, x4_binned_total = binning_data_more(x4, x3_binned_total, nbin, event_num,2)
	x4_patches_list, x5_binned_total = binning_data_more(x5, x4_binned_total, nbin, event_num,3)
	x5_patches_list=[]
	for bin in range(len(x5_binned_total)):
		x5_patches = binning_patches(x5_binned_total[bin],nbin)
		x5_patches_list.append(x5_patches)
		print(x5_patches)
	
	return x1_patches_list, x2_patches_list, x3_patches_list, x4_patches_list, x5_patches_list

### some repeated code used in find_c_t() ###
def binning_events(x, x_patches_list,num):
	mask_x_list     = []
	x_patches_order = []
	for i in range (num):
		mask_x_left  = (x >= x_patches_list[i][0]) & (x <  x_patches_list[i][1]) #find the index of the event in the a specific bin
		mask_x_right = (x >= x_patches_list[i][1]) & (x <= x_patches_list[i][2])
		mask_x_list.extend([mask_x_left, 
		                    mask_x_right])
		x_patches_order.extend([[x_patches_list[i][0],x_patches_list[i][1]],
		                        [x_patches_list[i][1],x_patches_list[i][2]]])
	return mask_x_list, x_patches_order

### find the C_T in the combination of each five bin of the phase-space variables ###
def find_c_t(bin_data,data,i_list,phase,type,event,event_n):
	x1_patches_list, x2_patches_list, x3_patches_list, x4_patches_list, x5_patches_list = find_bins(bin_data)
	data = np.transpose(data)
	x1,x2,x3,x4,x5,x6 = data[0],data[1],data[2],data[3],data[4],data[5]
	mask_x1_list, x1_patches_order = binning_events(x1, x1_patches_list, 1)
	mask_x2_list, x2_patches_order = binning_events(x2, x2_patches_list, 2)
	mask_x3_list, x3_patches_order = binning_events(x3, x3_patches_list, 4)
	mask_x4_list, x4_patches_order = binning_events(x4, x4_patches_list, 8)
	mask_x5_list, x5_patches_order = binning_events(x5, x5_patches_list, 16)
    #print the bin orientation
	c_p_list = []
	c_n_list = []
	i1,i2,i3,i4,i5 = i_list[0],i_list[1],i_list[2],i_list[3],i_list[4]
	for i in range (len(i1)):
	    with open('results_phase%.4f_%d/phase_space_regions_bin_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n),mode='a') as f:
	        f.write("(%.2f,%.2f)   (%.2f,%.2f)   (%.2f,%.2f)   (%.2f,%.2f)   (%.2f,%.2f)"%(
	                    x1_patches_order[int(i1[i])][0], x1_patches_order[int(i1[i])][1],
	                    x2_patches_order[int(i2[i])][0], x2_patches_order[int(i2[i])][1],
	                    x3_patches_order[int(i3[i])][0], x3_patches_order[int(i3[i])][1],
	                    x4_patches_order[int(i4[i])][0], x4_patches_order[int(i4[i])][1],
	                    x5_patches_order[int(i5[i])][0], x5_patches_order[int(i5[i])][1]) + "\n")
	    # find the events included in all the five bins
		mask_p = mask_x1_list[int(i1[i])] & mask_x2_list[int(i2[i])] & mask_x3_list[int(i3[i])] & mask_x4_list[int(i4[i])] & mask_x5_list[int(i5[i])] & (x6 > 0) #C_T>0
		mask_n = mask_x1_list[int(i1[i])] & mask_x2_list[int(i2[i])] & mask_x3_list[int(i3[i])] & mask_x4_list[int(i4[i])] & mask_x5_list[int(i5[i])] & (x6 < 0) #C_T<0
		# apply the event index to C_T
		c_p = x6[mask_p]
		c_n = x6[mask_n]
		c_p_list.append(c_p)
		c_n_list.append(c_n)
	return c_p_list,c_n_list

### find the A_T ###
def find_a_t(c_p_list,c_n_list):
	C_T_p = 1.0 * np.array([len(x) for x in c_p_list])
	C_T_n = 1.0 * np.array([len(x) for x in c_n_list])
	Asy   = (C_T_p-C_T_n) / (C_T_p+C_T_n)
	err_A = ((2 * C_T_n / (C_T_p + C_T_n) ** 2) ** 2 * C_T_p + (2 * C_T_p / (C_T_p + C_T_n) ** 2) ** 2 * C_T_n) ** 0.5
	return Asy,err_A
	
def line(x, b):
    return 0 * x +b

def main (program,type,event,phase,event_n):
	data      = np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n))
	data_conj = np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d_conj.txt'%(phase,event_n,type,event,event_n))

	#bin arrangement options
	i1 = [0,0,0,0,1,1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0,0,0,0,1,1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]   #paper
	i2 = [0,0,1,1,2,2, 3, 3, 0, 0, 1, 1, 2, 2, 3, 3, 0,0,1,1,2,2, 3, 3, 0, 0, 1, 1, 2, 2, 3, 3]
	i3 = [0,1,2,3,4,5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0,1,2,3,4,5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]
	i4 = [0,2,4,6,8,10,12,14,1, 3, 5, 7, 9, 11,13,15,0,2,4,6,8,10,12,14,1, 3, 5, 7, 9, 11,13,15]
	i5 = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,1,3,5,7,9,11,13,15,17,19,21,23,25,26,29,31]
 
	# i1 = [0,0,0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]    #Shyam
	# i2 = [0,0,0,0,0,0,0,0,1,1,1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]
	# i3 = [0,0,0,0,1,1,1,1,2,2,2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7]
	# i4 = [0,0,1,1,2,2,3,3,4,4,5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,10,11,11,12,12,13,13,14,14,15,15]
	# i5 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
	#  
	# i1 = [0,1,0,1,0,1,0,1,0,1,0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
	# i2 = [0,1,2,3,0,1,2,3,0,1,2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
	# i3 = [0,1,2,3,4,5,6,7,0,1,2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]
	# i4 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15]
	# i5 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
	#  
	# i1 = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]
	# i2 = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]
	# i3 = [0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7]
	# i4 = [0,2,4,6,8,10,12,14,1, 3, 5, 7, 9, 11,13,15,0,2,4,6,8,10,12,14,1, 3, 5, 7, 9, 11,13,15]
	# i5 = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,1,3,5,7,9,11,13,15,17,19,21,23,25,26,29,31]

	i_list = [i1,i2,i3,i4,i5]
	c_p_list,c_n_list = find_c_t(data,data,i_list,phase,type,event,event_n)
	c_p_list_conj,c_n_list_conj = find_c_t(data,data_conj,i_list,phase,type,event,event_n)
	# print(len(c_p_list))    #number of ps regions
	# print(len(c_n_list))
	# print([len(x) for x in c_p_list])   #number of C_T > 0 in each region
	# print([len(x) for x in c_n_list])   #number of C_T < 0 in each region
	# print([len(x) for x in c_p_list]+[len(x) for x in c_n_list])
	# print(np.sum([len(x) for x in c_p_list])) #total number of C_T > 0
	# print(np.sum([len(x) for x in c_n_list])) #total number of C_T < 0
	# print(np.sum([len(x) for x in c_p_list])+np.sum([len(x) for x in c_n_list])) #total number of C_T

	#calculation
	Asy,      err_A      = find_a_t (c_p_list,      c_n_list)
	Asy_conj, err_A_conj = find_a_t (c_p_list_conj, c_n_list_conj)
	a_cp    = 0.5 * (Asy - Asy_conj)
	err_acp = 0.5 * np.sqrt(err_A ** 2 + err_A_conj ** 2)

	for i in range (len(Asy)):
		with open('results_phase%.4f_%d/phase_space_regions_asy_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n),mode='a') as f:
			f.write("(%.2f +- %.2f)%%     (%.2f +- %.2f)%%     (%.2f +- %.2f)%%"%(
					 Asy[i]      * 100, err_A[i]      * 100, 
					 Asy_conj[i] * 100, err_A_conj[i] * 100,
					 a_cp[i]     * 100, err_acp[i]    * 100) + "\n")

	#fitting   
	ps_region  = np.arange(1,33)
	popt, pcov = curve_fit(line, ps_region, a_cp, sigma=err_acp)
	perr       = np.sqrt(np.diag(pcov))

	#print fit parameters and 1-sigma estimates
	print('fit parameter 1-sigma error')

	for i in range(len(popt)):
		print(str(popt[i])+'\pm'+str(perr[i]))

	# prepare confidence level curves
	nstd    = 5. # to draw 5sigma intervals
	popt_up = popt + nstd * perr
	popt_dw = popt - nstd * perr

	fit    = line(ps_region, *popt)
	fit_up = line(ps_region, *popt_up)
	fit_dw = line(ps_region, *popt_dw)

	#find chi2
	z       = (a_cp)/err_acp
	chi2    = np.sum(z ** 2)
	p_value = ROOT.TMath.Prob(chi2,32)

	z_fit       = (a_cp-fit)/err_acp
	chi2_fit    = np.sum(z_fit ** 2)
	p_value_fit = ROOT.TMath.Prob(chi2_fit,32)

	#plot
	fig, axs   = plt.subplots(1, 3, figsize=(15,5))
	(ax1, ax2, ax3) = axs
	fig.suptitle('Binned Analysis')

	ax1.errorbar(ps_region, Asy,      yerr=err_A,      fmt='ro' ,ms=4, capsize=2, lw=1)
	ax2.errorbar(ps_region, Asy_conj, yerr=err_A_conj, fmt='ro', ms=4, capsize=2, lw=1)
	ax3.errorbar(ps_region, a_cp,     yerr=err_acp,    fmt='ro', ms=4, capsize=2, lw=1)
	ax1.set(xlabel='Phase-Space Region', ylabel="$A_T$")
	ax2.set(xlabel='Phase-Space Region', ylabel="$\\bar{A}_T$")
	ax3.set(xlabel='Phase-Space Region', ylabel="$\\mathcal{A}_{\\mathcal{C}\\mathcal{P}}$")
	ax1.xaxis.set_label_coords(0.8, -0.11,  transform=ax1.transAxes)
	ax2.xaxis.set_label_coords(0.8, -0.11,  transform=ax2.transAxes)
	ax3.xaxis.set_label_coords(0.8, -0.11,  transform=ax3.transAxes)
	ax1.yaxis.set_label_coords(-0.12, 0.95, transform=ax1.transAxes)
	ax2.yaxis.set_label_coords(-0.12, 0.95, transform=ax2.transAxes)
	ax3.yaxis.set_label_coords(-0.12, 0.95, transform=ax3.transAxes)
	ax1.set_ylim([-0.42, 0.42])
	ax2.set_ylim([-0.42, 0.42])
	ax3.set_ylim([-0.11, 0.11])

	ax3.plot(ps_region, fit,         'b-', lw=2, label='best fit curve')
	ax3.plot(ps_region, np.zeros(32),'k--',lw=1, label='true curve')
	ax3.fill_between(ps_region, fit_up, fit_dw, alpha=0.25, label='$5\sigma$ interval')

	ax3.legend(loc='lower right', fontsize=12)
	ax3.text(0.01, 0.95, "No CPV: ${\chi}^2$=%.2f, p-value = %.1f%%"                     %(chi2,    p_value*100),     transform=ax3.transAxes)
	ax3.text(0.01, 0.88, "Fitted: $\;\;\;$${\chi}^2$=%.2f, p-value = %.1f%%"             %(chi2_fit,p_value_fit*100), transform=ax3.transAxes)
	ax3.text(0.2,  0.84, "($\\mathcal{A}_{\\mathcal{C}\\mathcal{P}}$ = %.4f $\pm$ %.4f)" %(popt[0], perr[0]),         transform=ax3.transAxes)

	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	fig.savefig("results_phase%.4f_%d/Binned_Analysis_factor%s%.1f_%d.png"%(phase,event_n,type,event,event_n))
    
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE        = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE       = float(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,PHASE,EVENT_N)
