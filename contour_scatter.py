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

# put all data into a list
yr = [
      #'0001','0002','0003','0004','0005',
      '0006','0007','0008','0009','0010'
     ]
mon = [
       '01','02','03','04','05','06','07','08','09','10','11','12'
      ]
Base = []
for i in range(len(yr)):
    for j in range(len(mon)):
        Base.append(xr.open_dataset
                     ('F:/Base.cam.h0.'+yr[i]+'-'+mon[j]+'.nc'))   
AClean = []
for i in range(len(yr)):
    for j in range(len(mon)):
        AClean.append(xr.open_dataset
                     ('F:/AClean.cam.h0.'+yr[i]+'-'+mon[j]+'.nc'))



# calculate concentration difference(variables indexed by m)
var = [
       'bc_a1_SRF','bc_a4_SRF','pom_a1_SRF','pom_a4_SRF','so4_a1_SRF',
       'so4_a2_SRF','so4_a3_SRF','soa_a1_SRF','soa_a2_SRF'
       #AEROD_v
      ]
nmon = len(yr)*len(mon)
ACleandata = [[]]*len(var)          
Basedata = [[]]*len(var)
conc_diff = np.zeros((9,192,288))
for m in range(len(var)):
    for n in range(len(Base)):
        AClean_data = pd.DataFrame(np.reshape(AClean[n][var[m]].data,(192,288))).fillna(0) *np.reshape(AClean[n]['PS'].data,(192,288))/np.reshape(AClean[n]['TS'].data,(192,288))/287.05*1e9 
        Base_data = pd.DataFrame(np.reshape(Base[n][var[m]].data,(192,288))).fillna(0) *np.reshape(Base[n]['PS'].data,(192,288))/np.reshape(Base[n]['TS'].data,(192,288))/287.05*1e9
        ACleandata[m].append(AClean_data)
        Basedata[m].append(Base_data)                 

ACleandata = np.array(ACleandata)
Basedata = np.array(Basedata)  

#############here!!!!!!        


 
lon = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lat.data      
stat = np.zeros((len(var),len(lat),len(lon)))
p = np.zeros((len(var),len(lat),len(lon)))

index_scatter = [[[],[]]]*len(var)
for a in range(len(var)):
    for b in range(len(lat)):
        for c in range(len(lon)):    
                stat[a,b,c],p[a,b,c] = stats.ttest_ind(
                Basedata[a,:,b,c],
                ACleandata[a,:,b,c],
                equal_var=False                       )
                if p[a,b,c]<=0.05:
                    index_scatter[a][1].append(lat[b])
                    index_scatter[a][0].append(lon[c]) 
# add 4 points to extend fig                    
    index_scatter[a][0].append(-180)
    index_scatter[a][1].append(90)
    index_scatter[a][0].append(-180)
    index_scatter[a][1].append(-90)
    index_scatter[a][0].append(180)
    index_scatter[a][1].append(90)
    index_scatter[a][0].append(180)
    index_scatter[a][1].append(-90)
                
    
    