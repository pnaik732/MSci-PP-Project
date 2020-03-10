import numpy as np
import matplotlib.pyplot as plt


def main(progname,phase,factor,type,event_n):
	data = np.transpose(np.loadtxt('results_phase%.4f_%d/asymmetries_%d.txt'%(phase,event_n,event_n)))

	event        = data[1]
	mask         = (event == factor)
	A_T          = data[3][mask]
	A_T_err      = data[4][mask]
	A_T_conj     = data[5][mask]
	A_T_err_conj = data[6][mask]
	a_cp         = data[7][mask]
	a_cp_err     = data[8][mask]

	f = plt.figure(figsize=(12,5))

	plt.subplot(1, 2, 1)
	plt.errorbar(event_n, A_T,yerr=A_T_err, fmt="bo",ms=4,label = "$A_T$ (Seed 7)",capsize=2,lw=1)
	plt.errorbar(event_n, A_T_conj,yerr=A_T_err_conj, fmt="go",ms=4,label = "$\\bar{A}_T$ (Seed 5)",capsize=2,lw=1)
	plt.errorbar(event_n, a_cp,yerr=err_a_cp, fmt="ro",ms=4,label = "$\mathcal{A}_{\mathcal{C}\mathcal{P}}$",capsize=2,lw=1)
	plt.title('Asymmetries with SPD Waves - All Resonances')
	plt.grid(True,linestyle='--',lw=0.5)
	plt.xlabel('Number of Events')
	plt.ylabel('Asymmetry')
	plt.ylim(-0.025,0.025)
	plt.legend()

	plt.subplot(1, 2, 2)

	plt.plot(event_n, A_T_err, ".-",c="b",label = "$A_T$",lw=1)
	plt.plot(event_n, A_T_err_conj, ".-",c="g",label = "$\\bar{A}_T$",lw=1)
	plt.plot(event_n, err_a_cp, ".-",c="r",label = "$\mathcal{A}_{\mathcal{C}\mathcal{P}}$",lw=1)
	plt.title('The Errors in Asymmetries')
	plt.grid(True,linestyle='--',lw=0.5)
	plt.xlabel('Number of Events')
	plt.ylabel('$\\Delta$(Asymmetry)')
	plt.legend()
	plt.tight_layout()

	plt.savefig('results_phase%.4f_%d/event_number_plots_factor%s%.1f_%d.png'%(phase,event_n,type,event,event_n))

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    PHASE       = float(sys.argv[1])
    FACTOR      = float(sys.argv[2])
    TYPE        = str(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME,PHASE,FACTOR,TYPE,EVENT_N)
    