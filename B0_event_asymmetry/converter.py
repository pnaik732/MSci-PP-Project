import uproot 
import pandas as pd
import sys
import os
def main(program,type,event,seed,phase,b):
    tree = uproot.open("phase%.4f/output_B0toDDbar0K+pi-_%s%.1f_s%d.root"%(phase,type,event,seed))['DalitzEventList']
    df0 = tree.pandas.df()
    if seed >= 0 and seed < int(b/10):
        df0.to_csv("output_phase%.4f/B0toDDbar0K+pi-_%s%.1f.txt"%(phase,type,event), sep=' ', mode='a')
        fin = open("output_phase%.4f/B0toDDbar0K+pi-_%s%.1f.txt"%(phase,type, event), "rt")
        data = fin.read()
        data = data.replace('~', 'p')
        data = data.replace('#', 'm')
        fin.close()
    
        fin = open("output_phase%.4f/B0toDDbar0K+pi-_%s%.1f.txt"%(phase,type,event), "wt")
        fin.write(data)
        fin.close()
    else:
        df0.to_csv("output_phase%.4f/B0toDDbar0K+pi-_%s%.1f.txt"%(phase,type,event), sep=' ', mode='a', header=False)

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE = str(sys.argv[1])
    EVENT        = float(sys.argv[2])
    SEED        = int(sys.argv[3])
    PHASE       = float(sys.argv[4])
    B           = int(sys.argv[4])
    main(PROGNAME, TYPE, EVENT,SEED,PHASE,B)
