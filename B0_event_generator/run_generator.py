#*****************************************************************************************
# Run the Generator under AmpGen/build/bin/B0_event_generator/
#
# B0toDDbar0K+pi-_spd_old.opt contains all the resonance with S,P and D waves
# edit B0toDDbar0K+pi-_spd_data.opt to specify the resonances want to use to generate events
#
# To change the relative P waves 
# specify a range in the factors of the amplitude
# specify a phase in radians
#
# written by Tailin Zhu
#*****************************************************************************************

import subprocess
import pandas as pd
import numpy as np

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  process.communicate()
  #stdout = process.communicate()[0]
  #print ('STDOUT:{}'.format(stdout))
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
event_type        = "spd1E" 

#for a range of amplitude factors and phases
# array1      = np.arange(-5,-1)      #ranges for the factors in the amplitudes of the p-wave resonance vertices
# array2      = np.arange(-1,1,0.2)
# array3      = np.arange(1,6)
# array_tot   = np.concatenate([array1,array2,array3])
# array_phase = np.arange(1.25*np.pi,2*np.pi,0.25*np.pi)    #ranges of the phase in the p-wave resonance vertices

#for amplitude factor=1 and phase=0:
array_tot     = np.arange(1,2) 
array_phase   = np.arange(0,1)

event_n = 10000        #total event number generated
event_n_generator = int(event_n/num_seed_in_a_set)
#=========================================================================================

df   = pd.read_csv('B0toDDbar0K+pi-_spd_data.opt',sep="\s+")
list = [1,4,7,10]  #index of the resonance in the data frame that want to make change with

for p in array_phase:
	run_cmd("mkdir generated_phase%.4f_%d"%(p,event_n))
	for i in array_tot:
		locals()['df_phase%.4f' % (p)] = np.array(df.phase)	
		locals()['df%.2f'       % (i)] = np.array(df.amp)
		for e in list:
			locals()['df_phase%.4f' % (p)][e] += p         #change the index specifying the phase varying resonance
			locals()['df%.2f'       % (i)][e] *= (10 ** i) #change the index specifying the amplitude varying resonance

		df2       = df.copy()
		df2.amp   = locals()['df%.2f'       % (i)]
		df2.phase = locals()['df_phase%.4f' % (p)]

		with open("B0toDDbar0K+pi-_spd_head.opt") as f:
			lines = f.readlines()
			with open("B0toDDbar0K+pi-_spd_new.opt", "w") as f1:
				f1.writelines(lines)
				f1.write("\n")
		df2.to_csv('B0toDDbar0K+pi-_spd_new.opt', sep=' ', mode='a', header=False, index=False, float_format='%.15f')

		for j in range (a,b,interval):
			print("calculating 1E%.1f; seed %d; phase %.4f"%(i,j,p))
			run_cmd("../Generator B0toDDbar0K+pi-_spd_new.opt --Seed=%d --nEvents=%d --Output=generated_phase%.4f_%d/output_B0toDDbar0K+pi-_factor%s%.1f_s%d_%d.root"     % (j,               event_n_generator, p, event_n, event_type, i, j, event_n_generator)) 
			run_cmd("../Generator B0toDDbar0K+pi-_spd_new.opt --Seed=%d --nEvents=%d --Output=generated_phase%.4f_%d/output_B0toDDbar0K+pi-_factor%s%.1f_s%d_%d_conj.root"% (j+conj_seed_add, event_n_generator, p, event_n, event_type, i, j+conj_seed_add, event_n_generator)) 

		a += 1
		b += 1
