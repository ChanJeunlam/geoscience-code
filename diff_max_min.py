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

# define var to be draw 



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
                     ('/public/home/gaojy/output/archive/TRA_AClean/atm/hist/TRA_AClean.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))


#read lon&lat
lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
# calc diff&data needed in t-test
nmon = len(yr)*len(mon)
aclean_tdata = []
base_tdata = []
diff = np.zeros((len(lat),len(lon)))
for o in range(nmon):
    aclean_ = np.reshape(aclean[o]['soa_a1_SRF'].data+aclean[o]['soa_a2_SRF'].data,(192,288))*np.reshape(aclean[o]['PS'].data,(192,288))/np.reshape(aclean[o]['TS'].data,(192,288))/287.05*1e9
    base_ = np.reshape(base[o]['soa_a1_SRF'].data+base[o]['soa_a2_SRF'].data,(192,288))*np.reshape(base[o]['PS'].data,(192,288))/np.reshape(base[o]['TS'].data,(192,288))/287.05*1e9
    #aclean_ = pd.DataFrame(np.reshape(aclean[o]['FSNT'].data,(192,288))-np.reshape(aclean[o]['FSNTC'].data,(192,288))-np.reshape(aclean[o]['FLNT'].data,(192,288))+np.reshape(aclean[o]['FLNTC'].data,(192,288))).fillna(0)
    #base_ = pd.DataFrame(np.reshape(base[o]['FSNT'].data,(192,288))-np.reshape(base[o]['FSNTC'].data,(192,288))-np.reshape(base[o]['FLNT'].data,(192,288))+np.reshape(base[o]['FLNTC'].data,(192,288))).fillna(0) 
    diff += (aclean_-base_)/nmon
diff = np.array(diff)
print(np.amax(diff))
print(np.amin(diff))