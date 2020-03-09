import numpy as np
import matplotlib.pyplot as plt

def asy_CP(A_T, err_A_T, A_Tbar, err_A_Tbar):
    a_CP = 0.5 * (A_T - A_Tbar)
    err_a_CP = np.sqrt((0.5 * err_A_T)**2+(0.5 * err_A_Tbar)**2)
    return a_CP, err_a_CP
    
def plt_fig(input,ylabel,ylim,phase,name,event_n):
    fig = input
	fig.suptitle('Psi(3770)K*(892)')
	for ax in fig1.get_axes():
		ax.set(xlabel='Factors', ylabel=ylabel)
		ax.grid(True,linestyle='--',lw=0.5)
		ax.set_xscale('log')
		ax.set_ylim([-ylim, ylim])
	
	for ax in fig1.get_axes():
		ax.label_outer()
	
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	fig.savefig('results_phase%.4f_%d/phase_plots_%s_%d.txt'%(phase,name,event_n))

def main(progname,phase,event_n):
	# array_phase1 = np.arange(-1.5*np.pi,0,0.25*np.pi)
	# array_phase2 = np.arange(1.75*np.pi,0,-0.25*np.pi)
	# array_phase  = np.concatenate([[0.0],array_phase1,array_phase2])

	v = 2
	h = int(len(phase)/v)
	fig1, axs1 = plt.subplots(v, h, figsize=(15,5),sharex='col', sharey='row',
							gridspec_kw={'hspace': 0, 'wspace': 0})
	(ax1, ax2, ax3, ax4, ax5, ax6, ax7), (ax8, ax9, ax10, ax11, ax12, ax13, ax14) = axs1


	fig2, axs2 = plt.subplots(v, h, figsize=(15,5),sharex='col', sharey='row',
							gridspec_kw={'hspace': 0, 'wspace': 0})
	(ax21, ax22, ax23, ax24, ax25, ax26, ax27), (ax28, ax29, ax30, ax31, ax32, ax33, ax34) = axs2

	fig3, axs3 = plt.subplots(v, h, figsize=(15,5),sharex='col', sharey='row',
							gridspec_kw={'hspace': 0, 'wspace': 0})
	(ax41, ax42, ax43, ax44, ax45, ax46, ax47), (ax48, ax49, ax50, ax51, ax52, ax53, ax54) = axs3


	i=1; j=21 ; k=41
	for p in phase:
		locals()['data_%.4f'        % (p)]  = np.transpose(np.loadtxt('results_phase%.4f_%d/asymmetries_%d.txt'%(p,event_n)))
		locals()['factor_%.4f'      % (p)]  = locals()['data_%.4f' % (p)][1]
		locals()['A_T_%.4f'         % (p)]  = locals()['data_%.4f' % (p)][3]
		locals()['A_T_err%.4f'      % (p)]  = locals()['data_%.4f' % (p)][4]
		locals()['A_T_conj%.4f'     % (p)]  = locals()['data_%.4f' % (p)][5]
		locals()['A_T_conj_err%.4f' % (p)]  = locals()['data_%.4f' % (p)][6]
		locals()['a_cp%.4f'         % (p)]  = locals()['data_%.4f' % (p)][7]
		locals()['a_cp_err%.4f'     % (p)]  = locals()['data_%.4f' % (p)][8]

		locals()['ax%d' % (i)].errorbar(locals()['factor_%.4f'      % (p)], 
										locals()['A_T_%.4f'         % (p)],
								   yerr=locals()['A_T_err%.4f'      % (p)],      fmt="ro",ms=4,capsize=2,lw=1)
		locals()['ax%d' % (i)].text(0.1, 0.9, "phase = %.2f$\pi$"   % (p/np.pi), transform=locals()['ax%d' % (i)].transAxes)
		locals()['ax%d' % (j)].errorbar(locals()['factor_%.4f'      % (p)], 
										locals()['A_T_conj%.4f'     % (p)],   
								   yerr=locals()['A_T_conj_err%.4f' % (p)],      fmt="ro",ms=4,capsize=2,lw=1)
		locals()['ax%d' % (j)].text(0.1, 0.9, "phase = %.2f$\pi$"   % (p/np.pi), transform=locals()['ax%d' % (i)].transAxes)
		locals()['ax%d' % (k)].errorbar(locals()['factor_%.4f'      % (p)], 
										locals()['a_cp%.4f'         % (p)], 
								   yerr=locals()['a_cp_err%.4f'     % (p)],      fmt="bo",ms=4,capsize=2,lw=1)
		locals()['ax%d' % (k)].text(0.1, 0.9, "phase = %.2f$\pi$"   % (p/np.pi), transform=locals()['ax%d' % (i)].transAxes)

		i+=1;j+=1;k+=1

	plt_fig(fig1,'$A_T$',0.043,phase,"A_T",event_n)
	plt_fig(fig2,'$\\bar{A}_T$',0.043,phase,"A_T_conj",event_n)
	plt_fig(fig3,'$\\mathcal{A}_{\\mathcal{C}\\mathcal{T}}$',0.043,phase,"a_cp",event_n)

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    PHASE       = float(sys.argv[1])
    EVENT_N     = int(sys.argv[1])
    main(PROGNAME,PHASE,EVENT_N)
    
