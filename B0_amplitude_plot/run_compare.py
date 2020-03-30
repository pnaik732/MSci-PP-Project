import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def run_cmd(cmd):
  process = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
  #process.communicate()
  stdout = process.communicate()[0]
  print ('STDOUT:{}'.format(stdout))
  process.stdout.close()
  
def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

decay_model = "B0toDDbar0K+pi-_spd.opt"
data_file   = "Event.root"
output_file = "Event_output_nondup_conj.root"
run_cmd("../SignalOnlyFitter %s --DataSample %s --Plots %s"%(decay_model,data_file,output_file))



number = ["01","02","03","12","13","23","012","013","023","123"]
for i in number:
	run_cmd("python compare.py %s %s"%(output_file,i))
 	
fig,axes = plt.subplots(3,4,figsize=(100,60))
axs = trim_axs(axes, 10)

for ax, no in zip(axs, number):
#     ax.set_title('markevery=%s' % str(case))
    ax.imshow(mpimg.imread("compare%s.png"%(no))); ax.set_axis_off()

fig.tight_layout()
fig.savefig("compare_all_conj.png")