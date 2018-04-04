from sys import argv
import numpy as np
import pandas as pd

script, input_matrix = argv

df1 = pd.read_csv(input_matrix, sep='\t', index_col=0)
print df1
