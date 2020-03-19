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
  #process.communicate()
  stdout = process.communicate()[0]
  print ('STDOUT:{}'.format(stdout))
  process.stdout.close()
  return " "

data0 = "data_old/Data_sig_tos_weights-Run1"; opt_cut0 = 0.9979

data1 = "data_new/Data_sig_tos_weights-Run1"; opt_cut1 = 0.9968
data2 = "data_new/Data_sig_tos_weights-Run2"; opt_cut2 = 0.9693
data3 = "data_new/Data_sig_tis_weights-Run1"; opt_cut3 = 0.9988
data4 = "data_new/Data_sig_tis_weights-Run2"; opt_cut4 = 0.9708

run_cmd("mkdir output")
run_cmd("mkdir output/data_new")
run_cmd("mkdir output/data_old")
run_cmd("python converter.py %s %.5f"         % (data0,opt_cut0));      
run_cmd("python converter.py %s %.5f"         % (data1,opt_cut1));     
run_cmd("python converter.py %s %.5f"         % (data2,opt_cut2));      
run_cmd("python converter.py %s %.5f"         % (data3,opt_cut3));      
run_cmd("python converter.py %s %.5f"         % (data4,opt_cut4));      print("Converting data file: done")  
# 
name0  = "run1_old"
run_cmd("mkdir results_%s"%(name0))
run_cmd("python CM_variables.py %s %s"        % (data0,name0)); print("calculate CM variables: done")
run_cmd("python plot_CM_variables.py %s"      % (name0));      print("plot CM variables: done")
run_cmd("python plot_fit_matplotlib.py %s"    % (name0));      print("fit invariant mass (matplotlib): done")
run_cmd("python plot_fit_ROOT.py %s"          % (name0));      print("fit invariant mass (ROOT): done")

name  = "run1ADDrun2_new"
run_cmd("mkdir results_%s"%(name))
run_cmd("python CM_variables.py %s %s"        % (data1,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data2,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data3,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data4,name)); print("calculate CM variables: done")

run_cmd("python plot_CM_variables.py %s"      % (name));      print("plot CM variables: done")
run_cmd("python plot_fit_matplotlib.py %s"    % (name));      print("fit invariant mass (matplotlib): done")
# run_cmd("python plot_fit_ROOT.py %s"          % (name));      print("fit invariant mass (ROOT): done")
