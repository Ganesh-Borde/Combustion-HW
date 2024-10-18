"""
Constant-pressure, adiabatic kinetics simulation.
"""

import sys
import numpy as np

import cantera as ct
import csv

gas = ct.Solution('gri30.xml')
fuel_species = 'CH4'
gas.set_equivalence_ratio(1.0, fuel_species, 'O2:1.0, N2:3.76')
gas.TP = 1600.0, ct.one_atm

r = ct.IdealGasConstPressureReactor(gas)

sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas, extra=['t'])
outfile = open('explosion_P_data.csv', 'w')
csvfile = csv.writer(outfile)

print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
npoints =100

for n in range(npoints):
    time += 1.e-5
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    print('%10.3e %10.f %10.f %14.6e' % (sim.time, r.T,
                                          r.thermo.P, r.thermo.u))
    XC=states.X[:,gas.species_index('O2')]
    csvfile.writerow([time, r.T, r.thermo.P, XC])
outfile.close()
print('Output written to file explosion_P_data.csv')
                                     
# Plot the results if matplotlib is installed.
# See http://matplotlib.org/ to get it.
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(states.t, states.T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(states.t, states.X[:,gas.species_index('OH')])
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(states.t, states.X[:,gas.species_index('H')])
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(states.t, states.X[:,gas.species_index('H2')])
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")
