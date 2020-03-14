import numpy as np
import sys

def find_a_t(c_p,c_n):
	C_T_p = 1.0 * len(c_p)
	C_T_n = 1.0 * len(c_n)
	Asy   = (C_T_p-C_T_n) / (C_T_p+C_T_n)
	err_A = ((2 * C_T_n / (C_T_p + C_T_n) ** 2) ** 2 * C_T_p + (2 * C_T_p / (C_T_p + C_T_n) ** 2) ** 2 * C_T_n) ** 0.5
	return Asy,err_A
	
def main(program, type,event,phase,event_n):
	CM_data      = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n)))
	CM_data_conj = np.transpose(np.loadtxt('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d_conj.txt'%(phase,event_n,type,event,event_n)))
	
	C_T          = CM_data[5];       C_T_conj     = CM_data_conj[5]
	C_T_p        = C_T[(C_T > 0)];   C_T_p_conj   = C_T_conj[(C_T_conj > 0)];
	C_T_n        = C_T[(C_T < 0)];   C_T_n_conj   = C_T_conj[(C_T_conj < 0)];
	
	A_T,      err_A_T      = find_a_t (C_T_p,      C_T_n)
	A_T_conj, err_A_T_conj = find_a_t (C_T_p_conj, C_T_n_conj)
	
	a_cp    = 0.5 * (A_T - A_T_conj)
	err_a_cp = 0.5 * np.sqrt(err_A_T ** 2 + err_A_T_conj ** 2)
	
	with open('results_phase%.4f_%d/asymmetries_%d.txt'%(phase,event_n,event_n),mode='a') as f:
		f.write("%.4f   %.1f    %d    %.5f   %.5f   %.5f   %.5f   %.5f   %.5f"%(
		        phase,  event,  event_n, A_T, err_A_T, A_T_conj, err_A_T_conj, a_cp, err_a_cp) + "\n")

	
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE        = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE       = float(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,PHASE,EVENT_N)
