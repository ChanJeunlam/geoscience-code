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
var = 'PRECC'
var_ = 'PRECL'



# put all data into a list
yr = [
      #'0001','0002','0003','0004','0005',
      '0006','0007','0008','0009','0010',
      '0011','0012','0013','0014','0015',
      '0016','0017','0018','0019','0020'
     ]
mon = [
       '01','02','12',
       #'03','04','05',
       #'06','07','08',
       #'09','10','11'
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
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))


#read lon&lat
lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
# calc diff&data needed in t-test
nmon = len(yr)*len(mon)
aclean_tdata = []
base_tdata = []
diff = np.zeros((len(lat),len(lon)))
for o in range(nmon):
    aclean_ = pd.DataFrame(np.reshape(aclean[o][var].data,(192,288))+np.reshape(aclean[o][var_].data,(192,288))).fillna(0)
    base_ = pd.DataFrame(np.reshape(base[o][var].data,(192,288))+np.reshape(base[o][var_].data,(192,288))).fillna(0) 
    aclean_tdata.append(aclean_)
    base_tdata.append(base_)
    diff += (aclean_-base_)/nmon  
aclean_tdata = np.array(aclean_tdata)
base_tdata = np.array(base_tdata)  
aclean_tdata_ = ((np.reshape(aclean_tdata,(15,3,192,288))).swapaxes(0,1)).mean(axis=0)
base_tdata_ = ((np.reshape(base_tdata,(15,3,192,288))).swapaxes(0,1)).mean(axis=0)  
# calc p_value in t-test&data needed in scatter      
stat = np.zeros((len(lat),len(lon)))
p_value = np.zeros((len(lat),len(lon)))
index_scatter = [[],[]]
for p in range(len(lat)):
    for q in range(len(lon)):    
            stat[p,q],p_value[p,q] = stats.ttest_ind(
            base_tdata_[:,p,q],
            aclean_tdata_[:,p,q],
            equal_var=False                   
            )
            if p_value[p,q]<=0.1:
                index_scatter[1].append(lat[p])
                index_scatter[0].append(lon[q])                  
index_scatter[0].append(-180)
index_scatter[1].append(90)
index_scatter[0].append(-180)
index_scatter[1].append(-90)
index_scatter[0].append(180)
index_scatter[1].append(90)
index_scatter[0].append(180)
index_scatter[1].append(-90)




# plot 
diff = np.array(diff)*60*60*24*1000
diff, lon = add_cyclic_point(diff, coord=lon)
# seperate colorbar   
#vmax = max(abs(np.reshape(diff,192*289)))
print(vmax)
#vmin = -vmax
norm = colors.Normalize(vmin=vmin, vmax=vmax)   
# add figure,coastline&borders
fig = plt.figure(figsize = (24,24))
proj = cartopy.crs.PlateCarree(central_longitude=0)   
ax1 = fig.add_subplot(1, 1, 1, projection=proj)
ax1.outline_patch.set_linewidth(3)
ax1.add_feature(cartopy.feature.COASTLINE,linewidth=2)
reader = Reader('/public/home/gaojy/plot/Data_ipynb/bou2_4p.shp')
provinces = cartopy.feature.ShapelyFeature(
    reader.geometries(), proj, edgecolor='k', facecolor='none')
ax1.add_feature(provinces, linewidth=2)
#ax1.add_feature(cartopy.feature.BORDERS,linewidth=4)



# set extent
#extent=[65,145,5,57]
#ax1.set_extent(extent)
# set lon&lat
ax1.set_xticks([-160,-120,-80,-40,0,40,80,120,160])
ax1.set_yticks([-90,-60,-30,0,30,60,90])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.tick_params(labelsize = 30)
# set title
ax1.set_title('precipitation-DJF',fontdict = {'fontsize' : 40})  
interval = [-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1]
# plot contour
h1 = plt.contourf(lon, lat, diff,
                  interval,
                  cmap = plt.cm.BrBG,
                  norm = norm,
                  extend = 'both')
# plot scatter
ax1.plot(index_scatter[0],index_scatter[1], markersize=1,
         marker='.',color='black',linestyle=' ',transform=proj)
# add colorbar    
cbar1 = fig.add_axes([0.15,0.23,0.72,0.02]) 
cb1 = plt.colorbar(h1, 
                   cax=cbar1,
                   ticks=interval,
                   orientation='horizontal')
cb1.ax.tick_params(labelsize=30)  
cb1.outline.set_linewidth(2)
cb1.set_label('mm day⁻¹',fontsize=28)
# save pic
plt.savefig('/public/home/gaojy/plot/diff-PREC-DJF.jpeg',
            #dpi=300, 
            bbox_inches = 'tight')

#plt.show()
