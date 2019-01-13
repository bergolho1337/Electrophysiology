#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

data = np.genfromtxt("PRd.txt",delimiter='\t',names=['t','V','ca'])

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
plot.show()

# Grafico de m
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Ca")
axis.set_xlabel('t')
axis.set_ylabel('Ca')
axis.plot(data['t'],data['ca'],c='g',label='ca')
leg = axis.legend()
plot.grid()
plot.ylim([0,0.0003])
plot.show()
