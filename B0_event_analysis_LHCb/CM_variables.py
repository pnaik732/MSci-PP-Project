import numpy as np
import uproot 
import pandas as pd
import ROOT
import sys

def boost_to_mother(particle1,particle2,particle3,particle4,motherparicle):
    boosttomother = -(motherparicle.BoostVector())
    particle1.Boost(boosttomother)
    particle2.Boost(boosttomother)
    particle3.Boost(boosttomother)
    particle4.Boost(boosttomother)

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
def main (program,data,name):
	df=pd.read_csv('output/%s.txt'%(data),sep="\s+")
	for i in range (len(df.index)):
		#extract database to form four-vectors (px,py,pz,E)
		locals()['B0_lv_%d'      % (i)] = ROOT.TLorentzVector(df.B0_PX[i]/1000,      df.B0_PY[i]/1000,      df.B0_PZ[i]/1000,      df.B0_PE[i]/1000)
		locals()['D0_lv_%d'      % (i)] = ROOT.TLorentzVector(df.D0_PX[i]/1000,      df.D0_PY[i]/1000,      df.D0_PZ[i]/1000,      df.D0_PE[i]/1000)
		locals()['Dbar0_lv_%d'   % (i)] = ROOT.TLorentzVector(df.D0bar_PX[i]/1000,   df.D0bar_PY[i]/1000,   df.D0bar_PZ[i]/1000,   df.D0bar_PE[i]/1000)
		locals()['Kplus_lv_%d'   % (i)] = ROOT.TLorentzVector(df.K_Kst0_PX[i]/1000,  df.K_Kst0_PY[i]/1000,  df.K_Kst0_PZ[i]/1000,  df.K_Kst0_PE[i]/1000)
		locals()['Piminus_lv_%d' % (i)] = ROOT.TLorentzVector(df.Pi_Kst0_PX[i]/1000, df.Pi_Kst0_PY[i]/1000, df.Pi_Kst0_PZ[i]/1000, df.Pi_Kst0_PE[i]/1000)

		locals()['sWeights_%d'   % (i)] = df.sWeights[i]
		locals()['ID_%d'         % (i)] = df.K_Kst0_ID[i]
		boost_to_mother(locals()['D0_lv_%d'% (i)],locals()['Dbar0_lv_%d' % (i)],locals()['Kplus_lv_%d' % (i)],locals()['Piminus_lv_%d' % (i)],locals()['B0_lv_%d' % (i)])


#       check if it boost correctly		
# 		print(locals()['B0_lv_%d' % (i)][0],locals()['B0_lv_%d' % (i)][1],locals()['B0_lv_%d' % (i)][2],locals()['B0_lv_%d' % (i)][3])
# 		locals()['B0_lv_%d' % (i)].Boost(-locals()['B0_lv_%d' % (i)].BoostVector())
# 		print(locals()['B0_lv_%d' % (i)][0],locals()['B0_lv_%d' % (i)][1],locals()['B0_lv_%d' % (i)][2],locals()['B0_lv_%d' % (i)][3])
# 
		##the sum of the four-momentum of the daughter particles and pairs
		locals()['momD0Dbar0_%d' % (i)]      = locals()['D0_lv_%d'    % (i)] + locals()['Dbar0_lv_%d'   % (i)]
		locals()['momKpPim_%d'   % (i)]      = locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
		locals()['momD0Dbar0KpPim_%d' % (i)] = locals()['D0_lv_%d' % (i)] + locals()['Dbar0_lv_%d' % (i)] + locals()['Kplus_lv_%d' % (i)] + locals()['Piminus_lv_%d' % (i)]
	
		#find the five CM variables
		locals()['invmass_D0Dbar0_%d'     % (i)] = locals()['momD0Dbar0_%d' % (i)].M()
		locals()['invmass_KpPim_%d'       % (i)] = locals()['momKpPim_%d'   % (i)].M()
		locals()['cosTheta_D0Dbar0_D0_%d' % (i)] = coshel(locals()['D0_lv_%d'    % (i)], locals()['momD0Dbar0_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
		locals()['cosTheta_KpPim_Kp_%d'   % (i)] = coshel(locals()['Kplus_lv_%d' % (i)], locals()['momKpPim_%d'   % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])
		locals()['Phi_%d'                 % (i)] = planeangle(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)], locals()['Piminus_lv_%d' % (i)], locals()['momD0Dbar0KpPim_%d' % (i)])

		#triple products
		
		#D0.(K+ X pi-)
		#locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['Kplus_lv_%d' % (i)], locals()['Piminus_lv_%d' % (i)], locals()['D0_lv_%d' % (i)])
		#K+.(D0 X Dbar0)
		locals()['C_Kpd_D0cDbar0_%d' % (i)] = tripleproduct(locals()['D0_lv_%d' % (i)], locals()['Dbar0_lv_%d' % (i)], locals()['Kplus_lv_%d' % (i)])


		#find the invariant mass for the four particles
		locals()['invmass_D0Dbar0KpPim_%d' % (i)] = locals()['momD0Dbar0KpPim_%d' % (i)].M()
		if locals()['ID_%d' % (i)] == 321:
			with open('results_%s/B0_CM_variables.txt'%(name),mode='a') as f:
				f.write("%.8f %.8f %.8f %.8f %.8f %.8f %.10f %.8f" % (locals()['Phi_%d'                  % (i)],
														              locals()['invmass_D0Dbar0_%d'      % (i)],
														              locals()['invmass_KpPim_%d'        % (i)],
														              locals()['cosTheta_D0Dbar0_D0_%d'  % (i)],
														              locals()['cosTheta_KpPim_Kp_%d'    % (i)],
														              locals()['C_Kpd_D0cDbar0_%d'       % (i)],
														              locals()['sWeights_%d'             % (i)],
														              locals()['invmass_D0Dbar0KpPim_%d' % (i)])
														              + "\n")
			
			with open("results_%s/B0toDDbarK_pi-_LHCb.txt"%(name),mode='a') as f2:
				f2.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f"%(df.D0_PE[i]/1000,      df.D0_PX[i]/1000,      df.D0_PY[i]/1000,      df.D0_PZ[i]/1000,      
			  																									 df.D0bar_PE[i]/1000,   df.D0bar_PX[i]/1000,   df.D0bar_PY[i]/1000,   df.D0bar_PZ[i]/1000,
																												 df.K_Kst0_PE[i]/1000,  df.K_Kst0_PX[i]/1000,  df.K_Kst0_PY[i]/1000,  df.K_Kst0_PZ[i]/1000,  
																												 df.Pi_Kst0_PE[i]/1000, df.Pi_Kst0_PX[i]/1000, df.Pi_Kst0_PY[i]/1000, df.Pi_Kst0_PZ[i]/1000, 
																												 df.sWeights[i])+ "\n")	
		if locals()['ID_%d' % (i)] == -321:	
			with open('results_%s/B0_CM_variables_conj.txt'%(name),mode='a') as f:
				f.write("%.8f %.8f %.8f %.8f %.8f %.8f %.10f %.8f" % (locals()['Phi_%d'                  % (i)],
														              locals()['invmass_D0Dbar0_%d'      % (i)],
														              locals()['invmass_KpPim_%d'        % (i)],
														              locals()['cosTheta_D0Dbar0_D0_%d'  % (i)],
														              locals()['cosTheta_KpPim_Kp_%d'    % (i)],
														             -locals()['C_Kpd_D0cDbar0_%d'       % (i)],
														              locals()['sWeights_%d'             % (i)],
														              locals()['invmass_D0Dbar0KpPim_%d' % (i)])
														              + "\n")
			
			with open("results_%s/B0toDDbarK_pi-_conj_LHCb.txt"%(name),mode='a') as f2:
				f2.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f"%(df.D0_PE[i]/1000,      df.D0_PX[i]/1000,      df.D0_PY[i]/1000,      df.D0_PZ[i]/1000,      
																													  df.D0bar_PE[i]/1000,   df.D0bar_PX[i]/1000,   df.D0bar_PY[i]/1000,   df.D0bar_PZ[i]/1000,
																													  df.K_Kst0_PE[i]/1000,  df.K_Kst0_PX[i]/1000,  df.K_Kst0_PY[i]/1000,  df.K_Kst0_PZ[i]/1000,  
																													  df.Pi_Kst0_PE[i]/1000, df.Pi_Kst0_PX[i]/1000, df.Pi_Kst0_PY[i]/1000, df.Pi_Kst0_PZ[i]/1000, 
																													  df.sWeights[i])+ "\n")
			
if __name__ == '__main__':
	PROGNAME    = sys.argv[0]
	DATA        = str(sys.argv[1])
	NAME        = str(sys.argv[2])
	main(PROGNAME, DATA,NAME)
