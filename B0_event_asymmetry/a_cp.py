import numpy as np
import matplotlib.pyplot as plt



def asy_CP(A_T, err_A_T, A_Tbar, err_A_Tbar):
    a_CP = 0.5 * (A_T - A_Tbar)
    err_a_CP = np.sqrt((0.5 * err_A_T)**2+(0.5 * err_A_Tbar)**2)
    return a_CP, err_a_CP


array_phase1 = np.arange(-1.5*np.pi,0,0.25*np.pi)
array_phase2 = np.arange(1.75*np.pi,0,-0.25*np.pi)
array_phase=np.concatenate([[0.0],array_phase1,array_phase2])

for p in array_phase:
    locals()['a_t%.4f' % (p)] = np.loadtxt("acp_p_phase%.4f.txt"%(p))
    locals()['a_t_%.4f' % (p)] = np.transpose(locals()['a_t_%.4f' % (p)])
    locals()['a_tbar%.4f' % (p)] = np.loadtxt("acp_p_phase%.4f_conj.txt"%(p))
    locals()['a_tbar%.4f' % (p)] = np.transpose(locals()['a_tbar%.4f' % (p)])
    locals()['a_cp%.4f' % (p)], locals()['err_a_cp%.4f' % (p)]= asy_CP(locals()['a_t%.4f' % (p)][1],
                                                                       locals()['a_t%.4f' % (p)][2],
                                                                       locals()['a_tbar%.4f' % (p)][1],
                                                                       locals()['a_tbar%.4f' % (p)][2])
    print("\n=======")
    print("phase = %.4f; amplitude factor = %.7f; \n A_T = %.6f pm %.6f; \n A_Tbar = %.6f pm %.6f; \n a_CP = %.6f pm %.6f" % (p,locals()['a_t%.4f' % (p)][0],locals()['a_t%.4f' % (p)][1],locals()['a_t%.4f' % (p)][2],locals()['a_tbar%.4f' % (p)][1],locals()['a_tbar%.4f' % (p)][2],locals()['a_cp%.4f' % (p)], locals()['err_a_cp%.4f' % (p)]))


plt.show()
    
