import uproot 
import pandas as pd
import sys

def main(program,type,event,seed,phase,num,event_n_gen,event_n):
    tree = uproot.open("generated_phase%.4f_%d/output_B0toDDbar0K+pi-_factor%s%.1f_s%d_%d.root"%(phase,event_n,type,event,seed,event_n_gen))['DalitzEventList']
    df0 = tree.pandas.df()
    if seed >= 0 and seed < num:
        df0.to_csv("output_phase%.4f_%d/B0toDDbar0K+pi-_factor%s%.1f_%d.txt"%(phase,event_n,type,event,event_n), sep=' ', mode='a')
        fin = open("output_phase%.4f_%d/B0toDDbar0K+pi-_factor%s%.1f_%d.txt"%(phase,event_n,type, event,event_n), "rt")
        data = fin.read()
        data = data.replace('~', 'p')
        data = data.replace('#', 'm')
        fin.close()
    
        fin = open("output_phase%.4f_%d/B0toDDbar0K+pi-_factor%s%.1f_%d.txt"%(phase,event_n,type,event,event_n), "wt")
        fin.write(data)
        fin.close()
    else:
        df0.to_csv("output_phase%.4f_%d/B0toDDbar0K+pi-_factor%s%.1f_%d.txt"%(phase,event_n,type,event,event_n), sep=' ', mode='a', header=False)

if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    TYPE = str(sys.argv[1])
    EVENT        = float(sys.argv[2])
    SEED        = int(sys.argv[3])
    PHASE       = float(sys.argv[4])
    B           = int(sys.argv[5])
    ENO_GEN     = int(sys.argv[6])
    ENO         = int(sys.argv[7]) 
    main(PROGNAME, TYPE, EVENT,SEED,PHASE,B,ENO_GEN,ENO)