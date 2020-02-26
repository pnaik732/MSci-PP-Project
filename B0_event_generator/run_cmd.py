import subprocess
import pandas as pd
import numpy as np

#B0toDDbar0K+pi-_spd.opt contains all the resonance
#edit B0toDDbar0K+pi-_p.opt to specify the resonances used to generate events

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  process.communicate()
  #stdout = process.communicate()[0]
  #print ('STDOUT:{}'.format(stdout))
  process.stdout.close()
  return " "

df = pd.read_csv('B0toDDbar0K+pi-_spd_data.opt',sep="\s+")


array1 = np.arange(-5,-1) #range of the factors
array2 = np.arange(-1,1,0.2)
array3 = np.arange(1,6)
array_tot = np.concatenate([array1,array2,array3])

array_phase=np.arange(-1.5*np.pi,2*np.pi,0.25*np.pi) #range of the phases

a = 0; b=3000 #the seed number (based on the event file used)
for p in array_phase:
    run_cmd("mkdir phase%.4f"%(p))
    for i in array_tot:
        locals()['df_phase%.4f' % (p)]     = np.array(df.phase)
        locals()['df_phase%.4f' % (p)][0] += p      #change the index specify the phase varying resonance
        locals()['df%.2f' % (i)]    = np.array(df.amp)
        locals()['df%.2f' % (i)][0] = locals()['df%.2f' % (i)][0] * (10 ** i) #change the index specify the amplitude varying resonance
    
        df2 = df.copy()
        df2.amp   = locals()['df%.2f' % (i)]
        df2.phase = locals()['df_phase%.4f' % (p)]
    
        with open("B0toDDbar0K+pi-_spd_head.opt") as f:
            lines = f.readlines()
            with open("B0toDDbar0K+pi-_p.opt", "w") as f1:
                f1.writelines(lines)
                f1.write("\n")
        df2.to_csv('B0toDDbar0K+pi-_p.opt', sep=' ', mode='a', header=False,index=False,float_format='%.15f')
        
        for j in range (a,b,int(b/10)):
            print("calculating 1E%.1f; seed %d; phase %.4f"%(i,j,p))
            run_cmd("../Generator B0toDDbar0K+pi-_p.opt --Seed=%d --nEvents=1000 --Output=phase%.4f/output_B0toDDbar0K+pi-_p1E%.1f_s%d.root"%(j,p,i,j)) 
            run_cmd("../Generator B0toDDbar0K+pi-_p.opt --Seed=%d --nEvents=1000 --Output=phase%.4f/output_B0toDDbar0K+pi-_p1E%.1f_s%d_conj.root"%(j+1000,p,i,j+1000)) 

        a += 1
        b += 1


