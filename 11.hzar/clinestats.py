import pandas as pd
from math import * 
import numpy as np
import matplotlib.pyplot as plt

def function(x,a,b):
    return 1/(1+exp(-(x-a)*4/b))

colnames = ["SNP","centre","centre_CI","width","width_CI","deltaF"]
data = pd.read_csv('pop_mafs/allsum.csv', header=1, names=colnames)
centre = data.centre.tolist()
width = data.width.tolist()
widthsort = data.width.tolist()
centresort = data.centre.tolist()

widthsort.sort()
centresort.sort()
lim = widthsort[int(0.05*len(width))]
print(lim)
print(widthsort[0])
print(np.mean(centre))
print(np.mean(width))

#5992.8637887433
#984.885533980606
#14952.193248940293
#16691.396224628497