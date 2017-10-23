# -*- coding: utf-8 -*-

# Paul's balloon calculator

# Based extensively on the work done by:
# Richard Meadows:  https://richardeoin.github.io/sp/
# Vincent E. Lally: https://opensky.ucar.edu/islandora/object/technotes%3A20

# With thanks to TomasTT7: https://github.com/TomasTT7/SuperpressureBalloonsNotebook/blob/master/Superpressure_Balloons.ipynb

import numpy
import matplotlib.pyplot as plt

def get_float(msg, default):
    ''' Ask the user for a float value using the provided message and default value '''
    try:
        prompt = 'Enter the ' + msg + ' (Default = ' + str(default) + '): '
        val = float(input(prompt))
    except:
        val = float(default)
    return val

def supertemp(alt):
    ''' Return the approximate supertemperature using a 3rd order fit to the typical values in Lally 1967 (Table 9) '''
    return ((-3e-12 * alt**3) + (1e-07 * alt**2) - (7e-05 * alt) + 7.0649)

print 'Balloon Calculator'
print

d = get_float('membrane density (g/m^2)',46.3) # Get the membrane density (g/m^2)
t = get_float('membrane thickness (microns)',44)/1000000 # Get the membrane thickness (microns), convert to m

# Ask the user for the balloon type:
# 'Mylar': https://en.wikipedia.org/wiki/Mylar_balloon_(geometry)
# Tubular 'Mylar': (approximation) assumes a cylindrical body between 'Mylar' end caps
# Rectangular ('Paper Bag'): https://en.wikipedia.org/wiki/Paper_bag_problem
balloon = raw_input('Enter the balloon type (M for \'Mylar\', T for Tubular \'Mylar\', R for rectangular) (Default = M): ')
if balloon != 'R' and balloon != 'M' and balloon != 'T': balloon = 'M'

if balloon == 'M': # 'Mylar' Balloon
    dia = get_float('balloon diameter including seams (m)',1.08) # Get the balloon diameter including seams
    seam = get_float('seam width (mm)', 12)/1000 # Get the heat weld seam width and convert to m
    extra = get_float('number of extra seams', 0) # Any extra seams - these are assumed to be the full diameter of the balloon
    # Calculate the mass of the empty balloon:
    # Density * ((double the area of the circular layers) + (double thickness of seam area))
    m = d * ((2 * numpy.pi * ((dia/2)**2)) + (extra * 2 * dia * seam))
    a = (dia / 2) - (seam) # Calculate the radius of the deflated balloon (radius - seam width)
    r = a / 1.311 # Approximate the radius of the inflated balloon
    v = 1000 * numpy.pi * a * (r**2) * 2 / 3 # Calculate the volume in L (not m^3)
    print 'Inflated balloon volume will be %.1fL' % v

elif balloon == 'T': # Tubular 'Mylar' Balloon
    dia = get_float('balloon width (diameter) including seams (m)',1.08) # Get the balloon diameter including seams
    length = get_float('total balloon length including seams (m)',1.08) # Get the balloon length including seams
    if length < dia: raise ValueError('Length must be >= Diameter!') # Check that length >= diameter
    seam = get_float('seam width (mm)', 12)/1000 # Get the heat weld seam width and convert to m
    extra = get_float('number of extra seams', 0) # Any extra seams - these are assumed to be the full length of the balloon
    # Calculate the mass of the empty balloon:
    # Density * ((double the area of the semi-circular ends plus the body) + (double thickness of seam area))
    m = d * ((2 * ((numpy.pi * ((dia/2)**2)) + ((length - dia) * dia))) + (extra * 2 * length * seam))
    a = (dia / 2) - (seam) # Calculate the radius of the deflated balloon (radius - seam width)
    r = a / 1.311 # Approximate the radius of the inflated balloon
    v = 1000 * numpy.pi * a * (r**2) * 2 / 3 # Calculate the volume in L of the 'Mylar' end caps
    # Add the volume of the body
    # Vmin assumes that the body cross section is a meridional hoop with area 2r^2
    # Hadzhilazova 2008 Equation 50
    # https://www.researchgate.net/publication/239784795_Once_More_the_Mylar_Balloon
    vmin = v + (1000 * (length - dia) * 2 * r**2)
    # Vmax assumes that the entire body is a cylinder with circumference 4a
    vmax = v + (1000 * (length - dia) * 4 * a**2 / numpy.pi)
    v = vmax # For the float calculations, assume the volume is vmax
    print 'Inflated balloon volume will be between %.1fL and %.1fL' %(vmin, vmax)

else: # Rectangular ('Paper Bag')
    w = get_float('balloon width including seams (m)',1.08) # Get the balloon width including seams
    h = get_float('balloon length including seams (m)',1.08) # Get the balloon length including seams
    seam = get_float('seam width (mm)', 12)/1000 # Get the heat weld seam width and convert to m
    extra = get_float('number of extra seams', 0) # Any extra seams - these are assumed to be the full length of the balloon
    # Calculate the mass of the empty balloon:
    # Density * ((double the area of the rectangle) + (double thickness of seam area))
    m = d * ((2 * w * h) + (extra * 2 * h * seam))
    w = w - (2 * seam) # Calculate the width of the inflatable balloon
    h = h - (2 * seam) # Calculate the length of the inflatable balloon
    v = 1000 * w**3 * ((h / (numpy.pi * w)) - (0.142 * (1 - 10**(-h/w)))) # Calculate the volume in L (not m^3)
    print 'Inflated balloon volume will be %.1fL' % v

print 'Balloon mass will be %.1fg' % m
pl = get_float('payload weight (g)',60) # Get the payload weight (g)
f = get_float('free lift (g)',20) # Get the free lift (g)
tl = m + pl + f # Calculate the required total lift on launch
print 'Total lift required at launch is %.1fg' % tl
print

# Ask the user for the gas type: hydrogen or helium
gas = raw_input('Enter the gas (H or He) (Default = He): ')
if gas != 'H' and gas != 'He': gas = 'He'
if gas == 'He':
    # Calculate the volume of gas required on launch to produce tl of lift
    gv = tl / (1.225 - 0.1786) # He (0.1786 g/L)  vs  Air (1.225 g/L)
    gm = gv * 0.1786 # Calculate the mass of the gas
else:
    # Calculate the volume of gas required on launch to produce tl of lift
    gv = tl / (1.225 - 0.08988) # H (0.08988 g/L)  vs  Air (1.225 g/L)
    gm = gv * 0.08988 # Calculate the mass of the gas

print
print 'Required volume of %s is %.1fL' % (gas,gv) # Print the required volume of gas
# Check that the balloon is large enough
if gv > v: raise ValueError('Required volume of %s exceeds balloon volume!' % gas)
print 'Balloon will be %.1f%% full at launch' % (100 * gv / v)
tm = m + pl + gm # Calculate the total mass: balloon + payload + gas
print 'Total balloon mass including gas is %.1fg' % tm

# Calculate the atmospheric temperature and pressure based on
# https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
altitude_11 = numpy.linspace(0,11000, num=11001) # 0-11km in 1m increments
temperature_11 = 15.04 - (0.00649 * altitude_11)
pressure_11 = 101.29 * numpy.power(((temperature_11 + 273.1) / 288.08),5.256)

altitude_25 = numpy.linspace(11001,25000, num=14000) # 11-25km in 1m increments
temperature_25 = numpy.ones(14000) * -56.46
pressure_25 = 22.65 * numpy.exp(1.73 - (0.000157 * altitude_25))

altitude = numpy.concatenate([altitude_11, altitude_25])
temperature = numpy.concatenate([temperature_11, temperature_25])
pressure = numpy.concatenate([pressure_11, pressure_25])
density = pressure / (0.2869 * (temperature + 273.1))

# PV = nRT
# P = pressure of the gas (Pascals)
# V = volume of the gas (m^3)
# n = amount of gas (mol)
# R = gas constant (8.3144598(48) J.mol^−1.K^−1)
# T = absolute temperature of the gas (K)

# Molar mass of He = 4.002602 g/mol
# Molar mass of H(2) = 1.00794 g/mol

if gas == 'He':
    mol = gm / 4.002602
else:
    mol = gm / 1.00794

print 'Balloon contains %.3f mol of %s' % (mol, gas)
print

# Calculate zero pressure altitude
for alt in range(len(altitude)):
    P = mol * 8.3144598 * (temperature[alt] + 273.16) / (v / 1000) # Calculate the pressure using P = nRT / V
    P = P / 1000 # Convert to kPa
    if P >= pressure[alt]: break
print 'Assuming the gas temperature is the same as that of the atmosphere:'
print 'Balloon will reach zero pressure differential at an altitude of %.0fm' % alt
print

# Calculate maximum altitude for zero balloon stretch (atmospheric_density * v == tm)
for den in range(len(density)):
    if density[den] <= tm / v: break

print 'Maximum altitude with zero stretch (Gamma = 1.0) is %.0fm' % altitude[den]
print 'Temperature at that altitude is %.1fC' % temperature[den]
print 'Pressure at that altitude is %.1fkPa' % pressure[den]
print

# Calculate the gas superpressure at maximum altitude for zero balloon stretch
P = mol * 8.3144598 * (temperature[den] + 273.16) / (v / 1000)
P = P / 1000 # Convert to kPa
print 'Assuming the gas temperature is the same as that of the atmosphere:'
print 'The gas pressure is %.1fkPa' % P
Pdiff = P - pressure[den]
# Show the superpressure in kPa and atmospheres (bar)
print 'Superpressure is %.1fkPa or %.3f atm' % (Pdiff, (Pdiff / 101.325))
print

# Calculate the typical supertemperature
stemp = supertemp(altitude[den])
P = mol * 8.3144598 * (temperature[den] + 273.16 + stemp) / (v / 1000)
P = P / 1000 # Convert to kPa
print 'The typical supertemperature at that altitude could be %.1fC' % stemp
print 'The gas pressure becomes %.1fkPa' % P
Pdiff = P - pressure[den]
print 'Superpressure becomes %.1fkPa or %.3f atm' % (Pdiff, (Pdiff / 101.325))
print

# Calculate maximum circumferential and longitudinal stress for the Mylar balloon
if balloon == 'M': # 'Mylar' Balloon
    cstress = (Pdiff * 1000 * r / t) / 1000000
    lstress = (Pdiff * 1000 * 0.599 * r / (2 * t)) / 1000000
    print 'Maximum circumferential stress is %.1fMPa' % cstress
    print 'Maximum longitudinal stress is %.1fMPa' % lstress

# Gamma is a measure of balloon stretch: volume / initial_volume
# For increasing gamma, calculate the new float altitude
gammas = [] # Gammas
alts = [] # Altitudes
sps = [] # Superpressures
for gamma in numpy.arange(1., 2.5, 0.01): # For gamma in the range 1 to 2.5
    gammas.append(gamma) # Append the gamma
    vol = v * gamma # Calculate the new volume
    for den in range(len(density)): # Find the new solution for atmospheric_density * vol == tm
        if density[den] <= tm / vol: break
    alts.append(altitude[den]) # Append the altitude
    # Calculate the gas pressure at that altitude taking supertemperature into account
    stemp = supertemp(altitude[den])
    P = mol * 8.3144598 * (temperature[den] + 273.16 + stemp) / (vol / 1000)
    P = P / 1000 # Convert to kPa
    Pdiff = P - pressure[den] # Calculate the superpressure
    sps.append(Pdiff) # Append the superpressure

# Plot altitude and superpressure vs gamma
f, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(gammas, alts)
ax2.plot(gammas, sps)
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.xlabel('Gamma (Volume / Initial Volume)')
ax1.set_title('Float Altitude (m) and Superpressure (kPa) vs Gamma')
ax1.grid(True)
ax2.grid(True)
plt.tight_layout()
plt.show()



   











    
