#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

data = np.genfromtxt("data.dat",delimiter=' ',names=['t','V','Pi','alpha','Nai','Ki'])

# Grafico de V
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Potencial transmembranico")
axis.set_xlabel('tempo (s)')
axis.set_ylabel('Potencial (mV)')
axis.plot(data['t'],data['V'],c='b',label='V')
leg = axis.legend(loc=4)
plot.grid()
plot.show()

# Grafico da pressao intracelular
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Pressao intracelular")
axis.set_xlabel('tempo (s)')
axis.set_ylabel('Pressao (atm)')
axis.plot(data['t'],data['Pi'],c='r',label='Pi')
leg = axis.legend(loc=4)
plot.grid()
plot.show()

# Grafico da relacao de volume intracelular
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Relacao volume intracelular")
axis.set_xlabel('tempo (s)')
axis.set_ylabel('alpha')
axis.plot(data['t'],data['alpha'],c='g',label='alpha')
leg = axis.legend(loc=4)
plot.grid()
plot.show()

# Grafico da concentracao de Na+ intracelular
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Concentracao intracelular Na+")
axis.set_xlabel('tempo (s)')
axis.set_ylabel('Concentracao (mmol)')
axis.plot(data['t'],data['Nai'],c='k',label='Nai')
leg = axis.legend(loc=2)
plot.grid()
plot.show()

# Grafico da concentracao de K+ intracelular
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("Concentracao intracelular K+")
axis.set_xlabel('tempo (s)')
axis.set_ylabel('Concentracao (mmol)')
axis.plot(data['t'],data['Ki'],c='k',label='Ki')
leg = axis.legend(loc=4)
plot.grid()
plot.show()

