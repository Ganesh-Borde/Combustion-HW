"""
Constant-pressure, adiabatic kinetics simulation.
"""

import sys
import csv
import numpy as np

import cantera as ct

gas = ct.Solution('gri30.cti')
air = ct.Solution('air.xml')

gas.TPX = 900.0, ct.one_atm, 'CH4:1,O2:2,N2:7.52'
npoints = 1000
r = ct.IdealGasReactor(gas)
env = ct.Reservoir(air)

# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = ct.Wall(r, env)
w.expansion_rate_coeff = 1.0e6  # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.area = 1.0

sim = ct.ReactorNet([r])
time = 0.0
times = np.zeros(npoints)
temps = np.zeros(npoints)
pressure = np.zeros(npoints)
data = np.zeros((gas.n_species,npoints))

# write output CSV file for importing into Excel
outfile = open('explosion.csv', 'w')
csvfile = csv.writer(outfile)
csvfile.writerow(['t (ms)','T (K)','P'] + gas.species_names)

print('%10s %10s %10s' % ('t [s]','T [K]','P [Pa]'))
for n in range(npoints):
    time += 1.e-2
    sim.advance(time)
    times[n] = time * 1e3  # time in ms
    temps[n] = r.thermo.T
    pressure[n] = r.thermo.P
    for i in range(gas.n_species):
        data[i,n] = r.thermo.X[i]
    csvfile.writerow([times[n], temps[n], pressure[n]] + list(data[:,n]))
    print('%10.3e %10.3f %10.3f' % (sim.time, r.T,
                                           r.thermo.P))

outfile.close()
print('Output written to file explosion.csv')
