import uproot 
import pandas as pd

root_file_name = "output_B0toDDbar0K+pi-.root"
txt_file_name = 'B0toDDbar0K+pi-.txt'
tree = uproot.open(file_name)['DalitzEventList']
df0 = tree.pandas.df()
df0.to_csv(txt_file_name, sep=' ', mode='a')
fin = open((txt_file_name, "rt")
data = fin.read()
data = data.replace('~', 'p')
data = data.replace('#', 'm')
fin.close()

fin = open((txt_file_name, "wt")
fin.write(data)
fin.close()
