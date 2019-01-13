#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

filename = "cell.dat"
print("Ploting file %s" % filename)
data = np.genfromtxt(filename,delimiter='\t',names=['t','V'])

# Grafico de V
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("V")
axis.set_xlabel('t')
axis.set_ylabel('V')
axis.plot(data['t'],data['V'],c='b',label='V')
leg = axis.legend()
plot.grid()
plot.ylim([-100,60])
plot.savefig(filename+"_v.pdf")
#plot.show()
