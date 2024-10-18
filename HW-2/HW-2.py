# %%


import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
#import csv

# Input variables
fuel_species = 'C3H8'
npoints = 1 # Number of points to consider
# Pressures of interest
P = [10, 100, 200, 300, 500, 1000, 3000, 1e4, 2e4, 3e4, 1e5, 5e5, 1e6, 1e7, 1e8]
# Temperatures explosion occurs for each pressure
Temps = np.zeros(len(P))
Temps[0] = 1000
exploded = False

i = -1

for Pi in P:
    i += 1
    for Ti in range(300,1500,10):
        
        #Display
        print('T %9.f P %11.2f' % (Ti, Pi))
        dt = 1 # Time step size [s]
        time = 0.0 # Initial time of interest [s]
        
        # Set up the gas
        gas = ct.Solution('gri30.yaml')
        gas.set_equivalence_ratio(1.0, fuel_species, 'O2:1.0, N2:3.76')
        gas.TP = Ti, Pi
        
        # Create reactor containing the gas
        reac = ct.IdealGasConstPressureReactor(gas)
        # Add the reactor to a new ReactorNet simulator
        sim = ct.ReactorNet([reac])
        
        states = ct.SolutionArray(gas)
        ##############################################################################
        # For loop to advance in time n times
        for n in range(npoints):
            time += dt
            # Advance sim to specified time
            sim.advance(time)
            states.append(reac.thermo.state)
            if np.abs(Ti-gas.T)>10:
                exploded = True
                break
            
        if (exploded):
            exploded = False
            Temps[i] = Ti
            break


plt.plot(Temps,P,'x-')
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.yscale('log')


# %%



fixed_temperature = 1100  # K
pressures = [20, 50,100]  # atm

# Initialize list to store induction times
induction_times = []


plt.figure(figsize=(10, 6))

for p in pressures:
    gas = ct.Solution('gri30.yaml')
    fuel_species = 'C3H8'
    gas.set_equivalence_ratio(1.0, fuel_species, 'O2:1.0, N2:3.76')
    
    # Set temperature and pressure
    gas.TP = fixed_temperature, p * ct.one_atm  # Convert atm to Pa
    
  
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])
    
    
    time = 0.0
    dt = 1e-05
    ntimes = 1000
    
    states = ct.SolutionArray(gas, extra=['t'])
    
    gradT = np.zeros(ntimes)
    T_before = fixed_temperature
    
    # Induction time tracking
    induction_time = None

    # Advance the simulation
    for n in range(ntimes):
        time += dt
        sim.advance(time)
        states.append(r.thermo.state, t=time * 1e3)  # in ms

        # Calculate the temperature gradient
        gradT[n] = np.abs(T_before - gas.T) / dt
        T_before = gas.T

    # Identify the induction time as the time when the temperature gradient is maximum
    induction_time = states.t[np.argmax(gradT)]

    
    plt.plot(np.log(states.t), states.T, label=f'Pressure={p} atm')

    
    induction_times.append(induction_time)

# Add plot labels and title for temperature evolution
plt.xlabel('log(Time) (ms)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Evolution Over Time at T=1000 K for Different Pressures')
plt.legend()
plt.grid()
plt.show()

# New figure for induction times
plt.figure(figsize=(10, 6))
plt.bar(pressures, induction_times, color='skyblue')
plt.xlabel('Pressure (atm)')
plt.ylabel('Induction Time (ms)')
plt.title('Induction Time at T=1000 K for Different Pressures')
plt.xticks(pressures)
plt.grid()
plt.show()


# %%
gas.TP = 1100.0, 20*ct.one_atm
r = ct.IdealGasConstPressureReactor(gas)

sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas, extra=['t'])


print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
npoints =1000

for n in range(npoints):
    time += 1.e-5
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    print('%10.3e %10.f %10.f %14.6e' % (sim.time, r.T,
                                          r.thermo.P, r.thermo.u))
    
# plt.plot(np.log(states.t), np.log(states.X[:,gas.species_index('OH')]),label=f'OH')


# plt.plot(np.log(states.t), np.log(states.X[:,gas.species_index('H')]),label=f'H')

# plt.plot(np.log(states.t), np.log(states.X[:,gas.species_index('C3H8')]),label=f'C3H8')

    

plt.plot(np.log(states.t), states.X[:,gas.species_index('OH')],label=f'OH')

plt.figure()
plt.plot(np.log(states.t), states.X[:,gas.species_index('H')],label=f'H')
plt.figure()
plt.plot(np.log(states.t), states.X[:,gas.species_index('C3H8')],label=f'C3H8')

plt.figure()
plt.plot(np.log(states.t), states.X[:,gas.species_index('CO')],label=f'CO')

plt.figure()
plt.plot(np.log(states.t), states.X[:,gas.species_index('NO')],label=f'NO')



# %%
gas.TP = 1100.0, 20*ct.one_atm  # Set temperature and pressure
r = ct.IdealGasConstPressureReactor(gas)

sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas, extra=['t'])

# Print header
print('%10s %10s %10s %14s' % ('t [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))

npoints = 1000

# Capture initial species concentrations for normalization
initial_OH = gas['OH'].X[0]
initial_H = gas['H'].X[0]
initial_C3H8 = gas['C3H8'].X[0]
initial_CO = gas['CO'].X[0]
initial_NO = gas['NO'].X[0]

# Reactor simulation
for n in range(npoints):
    time += 1.e-5
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    print('%10.3e %10.f %10.f %14.6e' % (sim.time, r.T, r.thermo.P, r.thermo.u))

# Normalized plots
plt.plot(np.log(states.t), states.X[:, gas.species_index('OH')] / initial_OH, label='OH')
plt.plot(np.log(states.t), states.X[:, gas.species_index('H')] / initial_H, label='H')
plt.plot(np.log(states.t), states.X[:, gas.species_index('C3H8')] / initial_C3H8, label='C3H8')
plt.plot(np.log(states.t), states.X[:, gas.species_index('CO')] / initial_CO, label='CO')
plt.plot(np.log(states.t), states.X[:, gas.species_index('NO')] / initial_NO, label='NO')

# Adding labels and legend
plt.xlabel('log(t)')
plt.ylabel('Normalized Mole Fraction')
plt.legend()
plt.show()

# %%



