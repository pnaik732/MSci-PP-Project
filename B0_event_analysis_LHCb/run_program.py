#*****************************************************************************************
# Run the all the analysis and plotting code under B0_event_analysis/
# converter.py           - convert the .pkl files to .txt files
# converter_root.py      - convert the .txt files to .root files
# CM_variables.py        - extract data and calculate CM variables and triple products
# plot_CM_variables.py   - plot CM variables histogram   
# plot_fit_matplotlib.py - plot fitted invariant masses for C_T>0 and C_T<0 and find
#                          find the asymmetries using method 1 and 2
# asymmetries.py         - find the asymmetries (method3)
# 
# @author: Tailin Zhu
#*****************************************************************************************

import subprocess
import numpy as np

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  #process.communicate()
  stdout = process.communicate()[0]
  print ('STDOUT:{}'.format(stdout))
  process.stdout.close()

data01 = "data_old/Data_sig_tos_weights-Run1"; opt_cut0 = 0.9979
data02 = "data_old/Data_sig_tos_weights-Run2"; opt_cut0 = 0.9979

data1 = "data_new/Data_sig_tos_weights-Run1"; opt_cut1 = 0.9968
data2 = "data_new/Data_sig_tos_weights-Run2"; opt_cut2 = 0.9693
data3 = "data_new/Data_sig_tis_weights-Run1"; opt_cut3 = 0.9988
data4 = "data_new/Data_sig_tis_weights-Run2"; opt_cut4 = 0.9708

#convert pkl file to txt file
run_cmd("mkdir output")

# run_cmd("mkdir output/data_old")
# run_cmd("python converter.py %s %.5f"         % (data01,opt_cut0));
# run_cmd("python converter.py %s %.5f"         % (data02,opt_cut0));  

run_cmd("mkdir output/data_new")     
run_cmd("python converter.py %s %.5f"         % (data1,opt_cut1));     
run_cmd("python converter.py %s %.5f"         % (data2,opt_cut2));      
run_cmd("python converter.py %s %.5f"         % (data3,opt_cut3));      
run_cmd("python converter.py %s %.5f"         % (data4,opt_cut4));      print("Converting data file: done")  

# make calculations
# name01  = "run1_old"
# run_cmd("mkdir results_%s"%(name01))
# run_cmd("python CM_variables.py %s %s"        % (data01,name01)); print("calculate CM variables: done")
# run_cmd("python plot_CM_variables.py %s"      % (name01));       print("plot CM variables: done")
# run_cmd("python plot_fit_matplotlib.py %s"    % (name01));       print("fit invariant mass (matplotlib): done")
# run_cmd("python asymmetries.py %s"            % (name01));       print("asymmetries from cut data: done")
# 
# name02 = "run2ADDrun2_old"
# run_cmd("mkdir results_%s"%(name02))
# run_cmd("python CM_variables.py %s %s"        % (data01,name02)); print("calculate CM variables: done")
# run_cmd("python CM_variables.py %s %s"        % (data02,name02)); print("calculate CM variables: done")
# run_cmd("python plot_CM_variables.py %s"      % (name02));       print("plot CM variables: done")
# run_cmd("python plot_fit_matplotlib.py %s"    % (name02));       print("fit invariant mass (matplotlib): done")
# run_cmd("python asymmetries.py %s"            % (name02));       print("asymmetries from cut data: done")

name  = "run1ADDrun2_new_nondup"
# name  = "run1ADDrun2_new_BW"
run_cmd("mkdir results_%s"%(name))
run_cmd("python CM_variables.py %s %s"        % (data1,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data2,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data3,name)); print("calculate CM variables: done")
run_cmd("python CM_variables.py %s %s"        % (data4,name)); print("calculate CM variables: done")

run_cmd("python plot_CM_variables.py %s"      % (name));      print("plot CM variables: done")
run_cmd("python plot_fit_matplotlib.py %s"    % (name));      print("fit invariant mass (matplotlib): done")
run_cmd("python asymmetries.py %s"            % (name));      print("asymmetries from cut data: done")
run_cmd("python converter_root.py %s"         % (name));      print("convert to .root file (amplitude analysis):done")