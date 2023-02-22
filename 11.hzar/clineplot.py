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

#hybrid index cline
#             Mean      SD
#center   15809.44 128.782
#width     8143.82 344.394
centre.append(15809.44)
width.append(8143.82)

widthsort.sort()
centresort.sort()
lim = widthsort[int(0.01*len(width))]
lim1 = centresort[int(0.1*len(centre))]
lim2 = centresort[int(0.9*len(centre))]
widthnew1 = []
widthnew2 = []
for i in centresort[0:int(0.1*len(centre))+1]:
    widthnew1.append(width[centre.index(i)])
widthnew1.sort()
lim3 = widthnew1[int(0.25*len(widthnew1))]
for i in centresort[int(0.9*len(centre))::]:
    widthnew2.append(width[centre.index(i)])
widthnew2.sort()
lim4 = widthnew2[int(0.25*len(widthnew2))]

xdata= np.linspace(0,30000,100)
for i in range(len(centre)):
    ydata = []
    for x in xdata:
        ydata.append(function(x,float(centre[i]),float(width[i])))
    if centre[i] == 15809.44 and width[i] == 8143.82:
        plt.plot(xdata,ydata,color="black",alpha=1,linewidth=1)
    elif width[i]<= lim:
        plt.plot(xdata, ydata, color='#808080', alpha=0.2, linewidth=0.3)
        #plt.plot(xdata, ydata, color='red', alpha=0.5, linewidth=0.8)
    #elif (centre[i] <= lim1 and width[i] <= lim3) or (centre[i] >= lim2 and width[i] <= lim4):
    #    plt.plot(xdata, ydata, 'g-',linewidth=0.7)
    else:
        plt.plot(xdata, ydata, color='#808080', alpha=0.2, linewidth=0.3)



plt.ylabel('allele frequency')
plt.xlabel('distance (m)')
plt.grid(True)
plt.savefig('all_withhi.pdf')  
