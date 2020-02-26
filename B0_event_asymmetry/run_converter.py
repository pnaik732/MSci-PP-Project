import subprocess
import os
import numpy as np
def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  process.communicate()
  process.stdout.close()
  return " "

a = 0; b=3000
array1 = np.arange(-5,-1)
array2 = np.arange(-1,1,0.2)
array3 = np.arange(1,6)
array_phase=np.arange(1.25*np.pi,2*np.pi,0.25*np.pi)
array_tot = np.concatenate([array1,array2,array3])
for p in array_phase:
    run_cmd("mkdir output_phase%.4f"%(p))
    for i in array_tot:
        for j in range (a,b,int(b/10)):
            print(p,j,i)
            run_cmd("python converter.py %s %.1f %d %.4f %d"%("p1E",i,j,p,b))
            #run_cmd("python converter_conj.py %s %.1f %d %.4f %d"%("p1E",i,j+1000,p,b))
        a +=1
        b +=1
    

    