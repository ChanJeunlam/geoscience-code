# import packages
import xarray as xr
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from cartopy.util import add_cyclic_point
import cartopy.mpl.ticker as mticker
import cartopy
import pandas as pd
import scipy.stats as stats
from cartopy.io.shapereader import Reader
import cmaps

# put all data into a list
yr = [
      #'0001','0002','0003','0004','0005',
      '0006','0007','0008','0009','0010',
      '0011','0012','0013','0014','0015',
      '0016','0017','0018','0019','0020'
     ]
mon = [
       '01','02','12',
       '03','04','05',
       '06','07','08',
       '09','10','11'
      ]
base = []
for m in range(len(yr)):
    for n in range(len(mon)):
        base.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))   
aclean = []
for m in range(len(yr)):
    for n in range(len(mon)):
        aclean.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/WST_AClean/atm/hist/WST_AClean.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))
#read lon&lat
lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
# calc diff&data needed in t-test
nmon = len(yr)*len(mon)
diff = np.zeros((len(lat),len(lon)))
for o in range(nmon):
    aclean_ = np.reshape(aclean[o]['soa_a1_SRF'].data+aclean[o]['soa_a2_SRF'].data,(192,288))*np.reshape(aclean[o]['PS'].data,(192,288))/np.reshape(aclean[o]['TS'].data,(192,288))/287.05*1e9
    base_ = np.reshape(base[o]['soa_a1_SRF'].data+base[o]['soa_a2_SRF'].data,(192,288))*np.reshape(base[o]['PS'].data,(192,288))/np.reshape(base[o]['TS'].data,(192,288))/287.05*1e9
    diff += (aclean_-base_)/nmon

flag_nc = xr.open_dataset('/public/home/gaojy/test/number_29_tag_regions.nc')
N = 0
conc_N = 0
for i in range(0,len(lat)):
    for j in range(0,len(lon)):
        if (flag_nc['region_code'].data[i,j]>=61 
            and flag_nc['region_code'].data[i,j]<=68):
            N += 1
            conc_N +=diff[i,j]
conc = conc_N/N

print(conc)
            
    
    
    
    
    
    
    
    
    
    
    