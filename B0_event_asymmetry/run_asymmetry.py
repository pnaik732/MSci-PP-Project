import subprocess
import numpy as np

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  process.communicate()
  process.stdout.close()
  return " "


array1 = np.arange(-5,-1)
array2 = np.arange(-1,1,0.2)
array3 = np.arange(1,6)
array_phase=np.arange(1.25*np.pi,2*np.pi,0.25*np.pi)
array_tot = np.concatenate([array1,array2,array3])

for p in array_phase:
    for i in array_tot:
        print(p,i)
        run_cmd("python asymmetry.py %s %.1f %.4f"%("p1E",i,p))
        #run_cmd("python asymmetry_conj.py %s %.1f %.4f"%("p1E",i,p))

