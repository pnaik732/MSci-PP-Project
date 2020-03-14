#*****************************************************************************************
# Run the all the analysis and plotting code under B0_event_analysis/
# converter.py           - convert ROOT data files to .csv files
# CM_variables.py        - extract data and calculate CM variables
# plot_CM_variables.py   - plot CM variables histogram   
# plot_fit_matplotlib.py - plot fitted invariant masses for C_T>0 and C_T<0
# plot_fit_ROOT.py       - plot fitted invariant masses for all events
# asymmetries.py         - find the triple product asymmetries and cp asymmetries
# binned_analysis.py     - find the asymmetries with binned data
# plot_phase.py          - plot with different phases
# plot_event.py          - plot with different event numbers
# 
# written by Tailin Zhu
#*****************************************************************************************

import subprocess
import os
import numpy as np

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  process.communicate()
  process.stdout.close()
  return " "

#=========================================================================================
#(copy this part to run_program.py in B0_event_analysis after any changes)

num_seed_in_a_set = 10    # number of seeds used to generate events for a event file  
a                 = 0
b                 = 10 * num_seed_in_a_set   # control the random seed (b-a >= num_seed_in_a_set)
conj_seed_add     = 1000  # add 1000 to seed number for generating the conjugate process 
interval          = int(b/num_seed_in_a_set)
interval_conj     = int(b/num_seed_in_a_set) + conj_seed_add 
event_type        = "spd_1E" 

#for a range of amplitude factors and phases
# array1      = np.arange(-5,-1)      #ranges for the factors in the amplitudes of the p-wave resonance vertices
# array2      = np.arange(-1,1,0.2)
# array3      = np.arange(1,6)
# array_tot   = np.concatenate([array1,array2,array3])
# array_phase = np.arange(1.25*np.pi,2*np.pi,0.25*np.pi)    #ranges of the phase in the p-wave resonance vertices

#for amplitude factor=1 and phase=0:
array_tot     = np.arange(0,1)  
array_phase   = np.arange(0,1)

event_n = 171300        #total event number generated
event_n_generator = int(event_n/num_seed_in_a_set)
#=========================================================================================

# for p in array_phase:
#     run_cmd("mkdir output_phase%.4f_%d"%(p,event_n))
#     for i in array_tot:
#         for j in range (a,b,interval):
#             print(p,j,i)
#             run_cmd("python converter.py %s %.1f %d %.4f %d %d %d"      % (event_type,i,j,              p,interval,     event_n_generator, event_n))
#             run_cmd("python converter_conj.py %s %.1f %d %.4f %d %d %d" % (event_type,i,j+conj_seed_add,p,interval_conj,event_n_generator, event_n))
#         a +=1
#         b +=1
#         
# print("Converting data file: done")  
 
f = 0.0  #choose a amplitude factor to plot asymmetries vs event numbers in plot_event.py
for p in array_phase:
	run_cmd("mkdir results_phase%.4f_%d"%(p,event_n))
	for i in array_tot:
		print(p,i)
		# run_cmd("python CM_variables.py %s %.1f %.4f %d"        % (event_type,i,p,event_n)); print("calculate CM variables: done")
# 		run_cmd("python CM_variables_conj.py %s %.1f %.4f %d"   % (event_type,i,p,event_n)); print("calculate CM variables_conj: done")
# 		run_cmd("python plot_CM_variables.py %s %.1f %.4f %d"   % (event_type,i,p,event_n)); print("plot CM variables: done")
		run_cmd("python plot_fit_matplotlib.py %s %.1f %.4f %d" % (event_type,i,p,event_n)); print("fit invariant mass (matplotlib): done")
		run_cmd("python plot_fit_ROOT.py %s %.1f %.4f %d"       % (event_type,i,p,event_n)); print("fit invariant mass (ROOT): done")
# 		run_cmd("python asymmetries.py %s %.1f %.4f %d"         % (event_type,i,p,event_n)); print("find asymmetries: done")
# 		run_cmd("python binned_analysis.py %s %.1f %.4f %d"     % (event_type,i,p,event_n)); print("binned analysis: done")
		#uncommand for multiple phases
		#run_cmd("python plot_phase.py %.4f %d"                  % (p,event_n));              print("phase varying plots: done")
		#uncommand for multiple event numbers
		#run_cmd("python plot_event.py %.4f %.1f %s %d"          % (p,f,event_type,event_n)); print("event number varying plots: done")
