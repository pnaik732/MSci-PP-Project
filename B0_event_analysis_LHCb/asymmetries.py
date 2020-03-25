import numpy as np
import sys

def find_a_t(c_p,c_n):
	C_T_p = 1.0 * len(c_p)
	C_T_n = 1.0 * len(c_n)
	Asy   = (C_T_p-C_T_n) / (C_T_p+C_T_n)
	err_A = ((2 * C_T_n / (C_T_p + C_T_n) ** 2) ** 2 * C_T_p + (2 * C_T_p / (C_T_p + C_T_n) ** 2) ** 2 * C_T_n) ** 0.5
	return Asy,err_A
	
def main(program, name):
	CM_data      = np.transpose(np.loadtxt('results_%s/B0_CM_variables.txt'%(name)))
	CM_data_conj = np.transpose(np.loadtxt('results_%s/B0_CM_variables_conj.txt'%(name)))
	
	mB0          = CM_data[7];       mB0_conj     = CM_data_conj[7]
	C_T          = CM_data[5];       C_T_conj     = CM_data_conj[5]
	C_T_p        = C_T[(C_T > 0) & (mB0 < 5.35) & (mB0 > 5.2)];   C_T_p_conj   = C_T_conj[(C_T_conj > 0) & (mB0_conj < 5.35) & (mB0_conj > 5.2)];
	C_T_n        = C_T[(C_T < 0) & (mB0 < 5.35) & (mB0 > 5.2)];   C_T_n_conj   = C_T_conj[(C_T_conj < 0) & (mB0_conj < 5.35) & (mB0_conj > 5.2)];
	
	A_T,      err_A_T      = find_a_t (C_T_p,      C_T_n)
	A_T_conj, err_A_T_conj = find_a_t (C_T_p_conj, C_T_n_conj)
	
	a_cp    = 0.5 * (A_T - A_T_conj)
	err_a_cp = 0.5 * np.sqrt(err_A_T ** 2 + err_A_T_conj ** 2)
	
	print("number of events (cut near the resonance): ",len(C_T_p)+len(C_T_n))
	with open('results_%s/invmass_fit_tp_parameters.txt'%(name),mode='a') as f:
		f.write("\n"+"=========================== Method 3 (from cuttings) =================================="+"\n")
		f.write("\n"+"Asymmetries (cut data)"+"\n")
		f.write("A_T       errA_T    A_T_conj errA_T_conj   a_cp     err_a_cp"+ "\n")
		f.write("%.4f   %.4f   %.4f   %.4f       %.4f   %.4f"%(
		         A_T,   err_A_T, A_T_conj, err_A_T_conj, a_cp, err_a_cp) + "\n")

	
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    NAME        = str(sys.argv[1])
    main(PROGNAME, NAME)
