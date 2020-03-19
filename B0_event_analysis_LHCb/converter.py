import pandas as pd
import sys

def main(progname, data, opt_cut):
	df = pd.read_pickle('../%s.pkl'%(data))
	df = df[df.NN_weights > opt_cut]
	df.to_csv('output/%s.txt'%(data), sep=' ', mode='a')
	print("total number of data before weight",len(df.index))
if __name__ == '__main__':
    PROGNAME    = sys.argv[0]
    DATA        = str(sys.argv[1])
    OPT_CUT     = float(sys.argv[2])
    main(PROGNAME, DATA,OPT_CUT)