import numpy as np
import pandas as pd
#from rootpy.vector import Vector2, Vector3, LorentzVector, Rotation, LorentzRotation
import sys
import ROOT

### find the triple product in the rest frame of mother particle###
def tripleproduct(particleA1,particleA2,particleB):
    particleA1v = particleA1.Vect()
    particleA2v = particleA2.Vect()
    particleBv  = particleB.Vect()
    cross       = particleA1v.Cross(particleA2v)
    dot         = particleBv.Dot(cross)
    return dot

### find the helicity angle in the rest frame  of daughter pairs###
def coshel(particle, parent, motherparticle):
    boosttoparent   = -(parent.BoostVector())
    particle.Boost(boosttoparent)
    motherparticle.Boost(boosttoparent)
    particle3       = particle.Vect()
    motherparticle3 = motherparticle.Vect()
    numerator       = particle3.Dot(motherparticle3)
    denominator     = (particle3.Mag()) * (motherparticle3.Mag())
    cosTheta        = numerator/denominator
    return cosTheta

### find the angle between the two planes in the rest frame of mother particle###
def planeangle(particleA1,particleA2,particleB1,particleB2,motherparticle):
    particleA1v = particleA1.Vect()
    particleA2v = particleA2.Vect()
    particleB1v = particleB1.Vect()
    particleB2v = particleB2.Vect()
    particleAv  = (particleA1 + particleA2).Vect()
    phi         = (particleA1v.Cross(particleA2v)).Angle(particleB1v.Cross(particleB2v))
    return phi
 
#========================== MAIN BODY
def main (program,type,event,phase,event_n):
	df=pd.read_csv("output_phase%.4f_%d/B0toDDbar0K+pi-_factor%s%.1f_%d.txt"%(phase,event_n,type,event,event_n),sep="\s+")
	for i in range (len(df.index)):
		#extract database to form four-vectors (px,py,pz,E)
		locals()['D0_lv_%d'      % (i)] = ROOT.TLorentzVector(df._1_D0_Px[i],    df._1_D0_Py[i],    df._1_D0_Pz[i],    df._1_D0_E[i])
		locals()['Dbar0_lv_%d'   % (i)] = ROOT.TLorentzVector(df._2_Dbar0_Px[i], df._2_Dbar0_Py[i], df._2_Dbar0_Pz[i], df._2_Dbar0_E[i])
		locals()['Kplus_lv_%d'   % (i)] = ROOT.TLorentzVector(df._3_Kp_Px[i],    df._3_Kp_Py[i],    df._3_Kp_Pz[i],    df._3_Kp_E[i])
		locals()['Piminus_lv_%d' % (i)] = ROOT.TLorentzVector(df._4_pim_Px[i],   df._4_pim_Py[i],   df._4_pim_Pz[i],   df._4_pim_E[i])
	
		#the sum of the four-momentum of the daughter-pair
		locals()['momD0Dbar0_%d' % (i)]      = locals()['D0_lv_%d'    % (i)] + locals()['Dbar0_lv_%d'   % (i)]
		locals()['momKpPim_%d'   % (i)]      = locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
		locals()['momD0Dbar0KpPim_%d' % (i)] = locals()['D0_lv_%d' % (i)] + locals()['Dbar0_lv_%d' % (i)] + locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
	
		#find the five CM variables
		locals()['invmass_D0Dbar0_%d'     % (i)] = locals()['momD0Dbar0_%d' % (i)].M()
		locals()['invmass_KpPim_%d'       % (i)] = locals()['momKpPim_%d'   % (i)].M()
		locals()['cosTheta_D0Dbar0_D0_%d' % (i)] = coshel(locals()['D0_lv_%d'    % (i)], locals()['momD0Dbar0_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
		locals()['cosTheta_KpPim_Kp_%d'   % (i)] = coshel(locals()['Kplus_lv_%d' % (i)], locals()['momKpPim_%d'   % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
		locals()['Phi_%d'                 % (i)] = planeangle(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)], locals()['Piminus_lv_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
	
		#find the triple products for K+ (D0 X Dbar0)
		locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)])
	
	#     #find the invariant mass for the four particles
	#     locals()['invmass_D0Dbar0KpPim_%d' % (i)] = locals()['momD0Dbar0KpPim_%d' % (i)].M()
	# 
	#     #check if the boost is correct
	#     locals()['momD0Dbar0_%d' % (i)].Boost(-locals()['momD0Dbar0_%d' % (i)].BoostVector())
	#     print(locals()['momD0Dbar0_%d' % (i)])
	
		#save results
		with open('results_phase%.4f_%d/B0_CM_variables_factor%s%.1f_%d.txt'%(phase,event_n,type,event,event_n),mode='a') as f:
			f.write("%.8f %.8f %.8f %.8f %.8f %.8f" % (locals()['Phi_%d'                 % (i)],
													   locals()['invmass_D0Dbar0_%d'     % (i)],
													   locals()['invmass_KpPim_%d'       % (i)],
													   locals()['cosTheta_D0Dbar0_D0_%d' % (i)],
													   locals()['cosTheta_KpPim_Kp_%d'   % (i)],
													   locals()['C_Kpd_D0cDbar0_%d'      % (i)])
													   + "\n")

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE        = str(sys.argv[1])
    EVENT       = float(sys.argv[2])
    PHASE       = float(sys.argv[3])
    EVENT_N     = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,PHASE,EVENT_N)
