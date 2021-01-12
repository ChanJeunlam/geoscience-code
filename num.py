import xarray as xr
from math import pi

emiss = xr.open_dataset('E:/2017_so4_a2_anthro-res_surface_mol.nc')

#emiss['AGR'].data = emiss['AGR'].data/6.022e23*115/1.77/(1/6*pi*(0.134e-4)**3)*6.022e26
#emiss['ENE'].data = emiss['ENE'].data/6.022e23*115/1.77/(1/6*pi*(0.261e-4)**3)*6.022e26
#emiss['IND'].data = emiss['IND'].data/6.022e23*115/1.77/(1/6*pi*(0.261e-4)**3)*6.022e26
#emiss['RCO'].data = emiss['RCO'].data/6.022e23*115/1.77/(1/6*pi*(0.0504e-4)**3)*6.022e26
#miss['SHP'].data = emiss['SHP'].data/6.022e23*115/1.77/(1/6*pi*(0.261e-4)**3)*6.022e26
#emiss['SLV'].data = emiss['SLV'].data/6.022e23*115/1.77/(1/6*pi*(0.134e-4)**3)*6.022e26
#emiss['TRA'].data = emiss['TRA'].data/6.022e23*115/1.77/(1/6*pi*(0.0504e-4)**3)*6.022e26
#emiss['WST'].data = emiss['WST'].data/6.022e23*115/1.77/(1/6*pi*(0.134e-4)**3)*6.022e26
#emiss['AGR'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['ENE'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['IND'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['RCO'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['SHP'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['SLV'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['TRA'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'
#emiss['WST'].attrs['units'] = '(particles/cm2/s)(molecules/mole)(g/kg)'

emiss.to_netcdf('E:/2017_num_so4_a2_anthro-res_surface_mol.nc','w','NETCDF4_CLASSIC')
