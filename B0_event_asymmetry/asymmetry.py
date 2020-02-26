import numpy as np
import ROOT 
from ROOT import gStyle,gPad
import uproot 
import pandas as pd
from rootpy.vector import Vector2, Vector3, LorentzVector, Rotation, LorentzRotation
import sys
#import matplotlib.pyplot as plt

    ### find the triple product in the rest frame of mother particle###
def tripleproduct(particleA1,particleA2,particleB):
    particleA1v = particleA1.Vect()
    particleA2v = particleA2.Vect()
    particleBv  = particleB.Vect()
    cross       = particleA1v.Cross(particleA2v)
    dot         = particleBv.Dot(cross)
    return dot
    
def main(program, type, event, phase):
    df=pd.read_table('output_phase%.4f/B0toDDbar0K+pi-_%s%.1f.txt'%(phase,type,event),sep="\s+")
    #========================== MAIN BODY

    C_T_positive_list = []
    C_T_negative_list = []

    for i in range (len(df.index)):
        #extract database to form four-vectors (px,py,pz,E)
        locals()['D0_lv_%d'      % (i)] = LorentzVector(df._1_D0_Px[i],    df._1_D0_Py[i],    df._1_D0_Pz[i],    df._1_D0_E[i])
        locals()['Dbar0_lv_%d'   % (i)] = LorentzVector(df._2_Dbar0_Px[i], df._2_Dbar0_Py[i], df._2_Dbar0_Pz[i], df._2_Dbar0_E[i])
        locals()['Kplus_lv_%d'   % (i)] = LorentzVector(df._3_Kp_Px[i],    df._3_Kp_Py[i],    df._3_Kp_Pz[i],    df._3_Kp_E[i])
        locals()['Piminus_lv_%d' % (i)] = LorentzVector(df._4_pim_Px[i],   df._4_pim_Py[i],   df._4_pim_Pz[i],   df._4_pim_E[i])
    
        locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)])

        if locals()['C_Kpd_D0cDbar0_%d' % (i)] > 0:
            C_T_positive_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])
        else:
            C_T_negative_list.append(locals()['C_Kpd_D0cDbar0_%d' % (i)])

    #Triple product asymmetry and the error
    C_T_p = 1.0 * len(C_T_positive_list)
    C_T_n = 1.0 * len(C_T_negative_list)
    Asy   = (C_T_p-C_T_n) / (C_T_p+C_T_n)
    err_A = ((2 * C_T_n / (C_T_p + C_T_n) ** 2) ** 2 * C_T_p + (2 * C_T_p / (C_T_p + C_T_n) ** 2) ** 2 * C_T_n) ** 0.5
    print(Asy, err_A)
    with open('acp_p_phase%.4f.txt'%(phase),mode='a') as f:
        f.write("%.7f %.6f %.6f"%(10**(event), Asy, err_A) + "\n")

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE       = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE      = float(sys.argv[3])
    main(PROGNAME, TYPE, EVENT, PHASE)

