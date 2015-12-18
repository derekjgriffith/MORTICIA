# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Tutorial in the Use of libRadtran

# <markdowncell>

# libRadtran is a collection of tools for computational radiative transfer in the Earth atmosphere. 
# This tutorial only shows the reading of uvspec input files and the reading of uvspec output files, as well as plotting of some outputs.
# 
# In order to run uvspec cases, a working installation of libRadtran is required, which is generally only possible on unix/linux machines.
# 
# For downloads and further information, go to http://www.libradtran.org

# <codecell>

# Use auto reload of librad for development purposes
%load_ext autoreload
%autoreload 1
%aimport librad
import librad
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# use latex for font rendering
mpl.rcParams['text.usetex'] = True  # Use TeX to format labels (takes a bit longer)
%matplotlib inline

# <codecell>

# Load a libRadtran example case
# Be default, any include files are expanded, creating a single set of option keywords
libRadCase=librad.Case(filename='examples/UVSPEC_AEROSOL.INP')

# <codecell>

# The printed version is identical to what is written in the libRadtran/uvspec input file
print libRadCase

# <codecell>

libRadCase.readout()  # Read uvspec output data for this case

# <codecell>

plt.figure(figsize=(10,8))
plt.plot(libRadCase.wvl, libRadCase.edir, libRadCase.wvl, libRadCase.edn, libRadCase.wvl, libRadCase.eup)
plt.xlabel('Wavlength [$nm$]')
plt.ylabel('Irradiance [$mW/m^2/sr/nm$]')
plt.legend(['Direct','Diffuse Downwelling','Total Upwelling'], loc='upper left')
plt.grid()

# <codecell>

# Thermal two-stream fast irradiances using the RODENTS solver
libRadCase=librad.Case(filename='examples/UVSPEC_RODENTS_ZOUT_THERMAL.INP')

# <codecell>

libRadCase.readout()

# <codecell>

plt.figure(figsize=(10,8))
plt.plot(libRadCase.zout, libRadCase.edn)
plt.title('RODENTS Solver : Downward Diffuse Irradiance')
plt.xlabel('Altitude Above Surface [km]')
plt.ylabel('Diffuse Irradiance [$mW/m^2/sr/nm$]')
plt.legend(['Fu Channel ' + str(chan) for chan in range(7,19)])
plt.xlim([0,30])
plt.grid()

# <codecell>

# Lidar solver example
# The lidar is at an altitude of 8.55 km, looking straight down
# Ground albedo is 0.2, causing sudden increase in returned signal at that range
# There is a cloud layer betweeen 2 and 4 km above ground
ss=librad.Case(filename='examples/UVSPEC_SSLIDAR.INP')

# <codecell>

ss.readout()  # Read the uvspec output

# <codecell>

plt.figure(figsize=(10,8))
plt.semilogy(ss.center_of_range, ss.number_of_photons)
plt.xlabel('Range from Laser [$km$]')
plt.ylabel('Returned Photons in Bin')
plt.grid()

# <codecell>

# This case shows use of the RPV (Rahman, Pinty, Verstraete [1993]) BRDF for a "plowed field"
# The wavelength is 400 nm (only)
x = librad.Case(filename='examples/UVSPEC_BRDF.INP')

# <codecell>

x.readout()

# <codecell>

x.uu.shape

# <codecell>

plt.figure(figsize=(10, 8))
plt.title('RPV BRDF for a Plowed Field')
plt.plot(x.umu, x.uu[:,:])  # Variable order is umu, phi, wvl, zout, nstokes
plt.xlabel('Cosine of Zenith Angle')
plt.ylabel('Transmittance at TOA')
plt.legend(['Azimuth ' + str(phi) + ' deg' for phi in x.phi])
plt.grid()

# <codecell>

# Here is a case with radiances in both umu and phi
y = librad.Case(filename='examples/UVSPEC_RADIANCES_ZOUT.INP')

# <codecell>

y.readout()

# <codecell>

print y.n_umu, y.n_phi, y.n_wvl, y.n_zout
print y.u0u.shape
y.u0u[:,1,0]

# <codecell>

# Plot radiances 
plt.figure(figsize=(10, 8), )
for iumu in range(y.n_umu):
    plt.plot(y.wvl, y.uu[iumu,0,:,0])  # order is umu, phi, wvl, zout, stokes
plt.xlabel('Wavelength [$nm$]')
plt.ylabel('Radiance [$mW/m^2/sr/nm$]')
plt.grid()

# <codecell>

plt.figure(figsize=(10, 8))
for iphi in range(y.n_phi):
    plt.plot(y.wvl, y.uu[2,iphi,:,1]) # order is umu, phi, wvl, zout, stokes
plt.xlabel('Wavelength [$nm$]')
plt.ylabel('Radiance [$mW/m^2/sr/nm$]')
plt.grid()

# <codecell>

# Solar irradiance profiles with old LOWTRAN model (fast, but not very accurate)
z = librad.Case(filename='examples/UVSPEC_PROFILES4.INP')

# <codecell>

z.readout()

# <codecell>

plt.figure(figsize=(10, 8))
plt.plot(z.wvl, z.edir.T)
plt.title('Direct Solar Irradiance (LOWTRAN Model)')
plt.legend(['Altitude ' + str(zout) + ' km' for zout in z.zout])
plt.xlabel('Wavelength [$nm$]')
plt.grid()

# <codecell>

# Two-stream irradiances using the fine spectral resolution REPTRAN model
lRep = librad.Case(filename='examples/UVSPEC_REPTRAN_SOLAR.INP')

# <codecell>

lRep.readout()

# <codecell>

plt.figure(figsize=(10,8))
plt.plot(lRep.wvl, lRep.edir, lRep.wvl, lRep.edn)
plt.title('REPTRAN Solar Band Model : Irradiance')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Irradiance [$mW/m^2/sr/nm$]')
plt.legend(['Direct Solar', 'Diffuse Downwelling'])
plt.grid()

# <codecell>

z = librad.Case(filename='examples/UVSPEC_RADIANCES_ZOUT.INP')

# <codecell>

z.readout()

# <codecell>

k = librad.Case(filename='examples/UVSPEC_REPTRAN_THERMAL.INP')

# <codecell>

k.readout()

# <codecell>

k.output_quantity

# <codecell>

c = librad.Case(filename='examples/UVSPEC_FILTER_SOLAR.INP')

# <codecell>

c

# <codecell>

c.readout()

# <codecell>

c.rad_units

# <codecell>

p = librad.Case(filename='examples/UVSPEC_AEROSOL_OPAC_POLaZOUT.INP')

# <codecell>

p.readout()

# <codecell>

b = librad.Case(filename='examples/UVSPEC_MC.INP')

# <codecell>

b.readout()

# <codecell>

b.fluxline

# <codecell>

b.edir, b.edn, b.eup, b.uavgdir, b.uavgdn, b.uavgup

# <codecell>

b.wvl

# <codecell>

# Reproduce the sunglint case plot for the Sentinel 3 OLCI band 13
sung = librad.Case(filename='examples/GUI/radiance_sentinel_sunglint/radiance_sentinel_sunglint.INP')

# <codecell>

sung.readout()

# <codecell>

sung

# <codecell>

plt.figure(figsize=(10,8))
plt.plot(np.rad2deg(np.arcsin(-sung.umu[::-1])), sung.uu[1::2],
         np.rad2deg(np.arcsin(sung.umu)), sung.uu[-2::-2])
plt.xlim([-90,90])
plt.ylim([0,110])
plt.title('Sunglint in Sentinel 3 OLCI Band 13')
plt.ylabel(r"Radiance [$mW/m^2/sr/nm$]", fontsize = 12)
plt.xlabel(r"Viewing nadir angle $\theta$ (degree)", fontsize = 12)
plt.grid()

# <codecell>

np.fliplr(-sung.umu)

# <codecell>

import numpy as np

# <codecell>

x = np.array([y,y,y,y,y])

# <codecell>

x.reshape((len(x),-1))

# <codecell>

y = np.vstack((x,x,x))

# <codecell>

dir(y)

# <codecell>

x.T * y.T

# <codecell>

a = np.array([y,y])

# <codecell>


# <codecell>


