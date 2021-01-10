# import packages
import scipy.stats as stats
import xarray as xr
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from cartopy.util import add_cyclic_point
import cartopy.mpl.ticker as mticker
import cartopy
import pandas as pd

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


var = [
       'bc_a1_SRF','bc_a4_SRF','pom_a1_SRF','pom_a4_SRF','so4_a1_SRF',
       'so4_a2_SRF','so4_a3_SRF','soa_a1_SRF','soa_a2_SRF',
       'AEROD_v'
       ]


AClean_vardata = [[],[],[],[],[],[],[],[],[],[]]           
Base_vardata = [[],[],[],[],[],[],[],[],[],[]]
for m in range(len(var)):
    for n in range(len(Base)):
        AClean_data = pd.DataFrame(np.reshape(AClean[n][var[m]].data,(192,288))).fillna(0) *np.reshape(AClean[n]['PS'].data,(192,288))/np.reshape(AClean[n]['TS'].data,(192,288))/287.05*1e9 
        Base_data = pd.DataFrame(np.reshape(Base[n][var[m]].data,(192,288))).fillna(0) *np.reshape(Base[n]['PS'].data,(192,288))/np.reshape(Base[n]['TS'].data,(192,288))/287.05*1e9
        AClean_vardata[m].append(np.array(AClean_data))
        Base_vardata[m].append(np.array(Base_data))                 

AClean_vardata = np.array(AClean_vardata)
Base_vardata = np.array(Base_vardata)       
stat = np.zeros((10,192,288))
p = np.zeros((10,192,288))
lon = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lat.data
index_scatter = [[[],[]]]*10
for a in range(10):
    for b in range(192):
        for c in range(288):    
                stat[a,b,c],p[a,b,c] = stats.ttest_ind(
                Base_vardata[a,:,b,c],
                AClean_vardata[a,:,b,c],
                equal_var=False                       )
                if p[a,b,c]<=0.05:
                    index_scatter[a][0].append(lat[b])
                    index_scatter[a][1].append(lon[c])
                    
index_scatter[0][0].append(-180)
index_scatter[0][1].append(90)
index_scatter[0][0].append(-180)
index_scatter[0][1].append(-90)
index_scatter[0][0].append(180)
index_scatter[0][1].append(90)
index_scatter[0][0].append(180)
index_scatter[0][1].append(-90)
                
# add figure,coastline&borders
fig = plt.figure(figsize = (8,6))
proj = cartopy.crs.PlateCarree(central_longitude=0)   
ax1 = fig.add_subplot(1, 1, 1, projection=proj)
ax1.add_feature(cartopy.feature.COASTLINE)
ax1.add_feature(cartopy.feature.BORDERS)

# set extent
# extent=[-180,180,-90,80]
# ax1.set_extent(extent)

# set lon&lat
ax1.set_xticks([-160,-120,-80,-40,0,40,80,120,160])
ax1.set_yticks([-90,-60,-30,0,30,60,90])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.set_title(var[m],fontdict = {'fontsize' : 18})

#ax1.plot(index_scatter[0][0],index_scatter[0][1], markersize=5,marker='o',linestyle='',color='#3b3b3b',transform=proj)
ax1.scatter(index_scatter[0][0], index_scatter[0][1],s = 2,color = 'black'
            #markersize=0.1,marker='o',linestyle='',color='black'
            )

plt.show()
