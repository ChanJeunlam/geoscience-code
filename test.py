import xarray as xr
import numpy as np
from math import pi


nc2 = xr.open_dataset('E:/2013_so4_a1_anthro-ene_vertical_mol.nc')
nc4 = xr.open_dataset('E:/2013_num_so4_a1_anthro-ene_vertical_mol.nc')



emiss2 = np.sum(np.sum(np.sum(np.sum(nc2.ENE.data+nc2.IND.data
          ,0),0),0),0)/6.022e23*115/1.77/(1/6*pi*(0.261e-4)**3)*6.022e26
emiss4 = np.sum(np.sum(np.sum(np.sum(nc4.num_so4_a1_ene_ind.data,0),0),0),0)
