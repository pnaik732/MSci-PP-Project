import pandas as pd
import numpy as np
import uproot
import sys

def convert(file,df):
	file.newtree("DalitzEventList",{'_1_D0_E':np.float64, 
	                              '_1_D0_Px':np.float64, 
	                              '_1_D0_Py':np.float64,
	                              '_1_D0_Pz':np.float64,
	                              '_2_Dbar0_E':np.float64, 
	                              '_2_Dbar0_Px':np.float64,
	                              '_2_Dbar0_Py':np.float64,
	                              '_2_Dbar0_Pz':np.float64,
	                              '_3_K~_E':np.float64,
	                              '_3_K~_Px':np.float64,
	                              '_3_K~_Py':np.float64,
	                              '_3_K~_Pz':np.float64,
	                              '_4_pi#_E':np.float64,
	                              '_4_pi#_Px':np.float64,
	                              '_4_pi#_Py':np.float64,
	                              '_4_pi#_Pz':np.float64,
	                              'weight':np.float64},"DalitzEventList")
	file["DalitzEventList"].extend({'_1_D0_E':df[0], 
	                              '_1_D0_Px':df[1]  , 
	                              '_1_D0_Py':df[2] ,
	                              '_1_D0_Pz':df[3]  ,
	                              '_2_Dbar0_E':df[4] , 
	                              '_2_Dbar0_Px':df[5]  ,
	                              '_2_Dbar0_Py':df[6]  ,
	                              '_2_Dbar0_Pz':df[7]  ,
	                              '_3_K~_E':df[8] ,
	                              '_3_K~_Px':df[9] ,
	                              '_3_K~_Py':df[10] ,
	                              '_3_K~_Pz':df[11]  ,
	                              '_4_pi#_E':df[12]  ,
	                              '_4_pi#_Px':df[13] ,
	                              '_4_pi#_Py':df[14] ,
	                              '_4_pi#_Pz':df[15] ,
	                              'weight':df[16]} )


def main(progname, name):
	df     =np.transpose(np.loadtxt("results_%s/B0toDDbarK_pi-_LHCb.txt"     %(name)))
	df_conj=np.transpose(np.loadtxt("results_%s/B0toDDbarK_pi-_conj_LHCb.txt"%(name)))

# 	df.columns = ['_1_D0_E', '_1_D0_Px', '_1_D0_Py', '_1_D0_Pz', '_2_Dbar0_E', '_2_Dbar0_Px', '_2_Dbar0_Py','_2_Dbar0_Pz', '_3_K~_E', '_3_K~_Px', '_3_K~_Py', '_3_K~_Pz', '_4_pi#_E', '_4_pi#_Px', '_4_pi#_Py', '_4_pi#_Pz', 'weight']
# 	df_conj.columns = ['_1_Dbar0_E', '_1_D0bar_Px', '_1_Dbar0_Py', '_1_Dbar0_Pz', '_2_D0_E', '_2_D0_Px', '_2_D0_Py','_2_D0_Pz', '_3_K#_E', '_3_K#_Px', '_3_K#_Py', '_3_K#_Pz', '_4_pi~_E', '_4_pi~_Px', '_4_pi~_Py', '_4_pi~_Pz', 'weight']
# 	df.insert(0, "entry", np.arange(0,len(df.index)), True)
# 	df_conj.insert(0, "entry", np.arange(0,len(df_conj.index)), True)
# 
# 	df.to_csv("results_%s/B0toDDbarK_pi-_LHCb.txt"     %(name), sep=' ', mode='w',index=False)
# 	df_conj.to_csv("results_%s/B0toDDbarK_pi-_conj_LHCb.txt"     %(name), sep=' ',mode='w',index=False)
# 	
	file = uproot.recreate("results_%s/Event.root"%(name))
	file_conj = uproot.recreate("results_%s/Event_conj.root"%(name))
	convert(file,df)
	convert(file_conj,df_conj)
	
	test = uproot.open("results_%s/Event.root"%(name))['DalitzEventList']
	df0 = test.pandas.df()
	df0.to_csv("results_%s/test.txt"%(name), sep=' ', mode='a')
	#df = df.rename_axis('entry', axis=1)
	#df_conj = df_conj.rename_axis('entry', axis=1)

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    NAME        = str(sys.argv[1])
    main(PROGNAME, NAME)