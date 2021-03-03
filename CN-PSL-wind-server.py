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
var = 'PSL'  
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
                     ('/public/home/gaojy/output/archive/AClean/atm/hist//AClean.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))   
aclean = []
for m in range(len(yr)):
    for n in range(len(mon)):
        aclean.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean_O3/atm/hist//AClean_O3.cam.h0.'+yr[m]+'-'+mon[n]+'.nc'))



lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist//AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist//AClean.cam.h0.0001-01.nc')).lat.data
lon_uv = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist//AClean.cam.h0.0001-01.nc')).lon.data
lat_uv = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist//AClean.cam.h0.0001-01.nc')).lat.data
# calc diff&data needed in t-test
nmon = len(yr)*len(mon)
aclean_tdata = []
base_tdata = []
aclean_u_tdata = []
base_u_tdata = []
aclean_v_tdata = []
base_v_tdata = []
diff = np.zeros((len(lat),len(lon)))
diff_u = np.zeros((len(lat),len(lon)))
diff_v = np.zeros((len(lat),len(lon)))
for o in range(nmon):
    aclean_ = pd.DataFrame(np.reshape(aclean[o][var].data,(192,288))).fillna(0)
    base_ = pd.DataFrame(np.reshape(base[o][var].data,(192,288))).fillna(0) 
    aclean_tdata.append(aclean_)
    base_tdata.append(base_)
    diff += (aclean_-base_)/nmon  
    aclean_u = np.reshape(aclean[o]['U'][0,25].data,(192,288))
    base_u = np.reshape(base[o]['U'][0,25].data,(192,288))
    aclean_u_tdata.append(aclean_u)
    base_u_tdata.append(base_u)
    diff_u += (aclean_u-base_u)/nmon 
    aclean_v = np.reshape(aclean[o]['V'][0,25].data,(192,288))
    base_v = np.reshape(base[o]['V'][0,25].data,(192,288))
    aclean_v_tdata.append(aclean_v)
    base_v_tdata.append(base_v)
    diff_v += (aclean_v-base_v)/nmon
aclean_tdata_mean = ((np.reshape(np.array(aclean_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0)
base_tdata_mean = ((np.reshape(np.array(base_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0) 
aclean_u_tdata_mean = ((np.reshape(np.array(aclean_u_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0)
base_u_tdata_mean = ((np.reshape(np.array(base_u_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0) 
aclean_v_tdata_mean = ((np.reshape(np.array(aclean_v_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0)
base_v_tdata_mean = ((np.reshape(np.array(base_v_tdata),(15,12,192,288))).swapaxes(0,1)).mean(axis=0) 



# calc p_value in t-test&data needed in scatter      
stat = np.zeros((len(lat),len(lon)))
p_value = np.zeros((len(lat),len(lon)))
stat_u = np.zeros((len(lat),len(lon)))
p_value_u = np.zeros((len(lat),len(lon)))
stat_v = np.zeros((len(lat),len(lon)))
p_value_v = np.zeros((len(lat),len(lon)))
index_x = []
index_y = []
index_uv_x = []
index_uv_y = []
index_uv_u = []
index_uv_v = []
for p in range(len(lat)):
    for q in range(len(lon)):    
            stat[p,q],p_value[p,q] = stats.ttest_ind(
            base_tdata_mean[:,p,q],
            aclean_tdata_mean[:,p,q],
            equal_var=False                   
            )
            stat_u[p,q],p_value_u[p,q] = stats.ttest_ind(
            base_u_tdata_mean[:,p,q],
            aclean_u_tdata_mean[:,p,q],
            equal_var=False                   
            )
            stat_v[p,q],p_value_v[p,q] = stats.ttest_ind(
            base_v_tdata_mean[:,p,q],
            aclean_v_tdata_mean[:,p,q],
            equal_var=False                   
            )
            if p_value[p,q]<=0.1:
                index_y.append(lat[p])
                index_x.append(lon[q])   
            #if p_value_u[p,q]>0.1 and p_value_v[p,q]>0.1:
            #    diff_u[p,q] = 0
            #   diff_v[p,q] = 0
                #index_uv_y.append(lat[p])
                #index_uv_x.append(lon[q])
                #index_uv_u.append(diff_u[p,q])
                #index_uv_v.append(diff_v[p,q])
index_x.append(-180)
index_y.append(90)
index_x.append(-180)
index_y.append(-90)
index_x.append(180)
index_y.append(90)
index_x.append(180)
index_y.append(-90)



# set contour interval
#interval = [-10.5,-7.5,-4.5,-1.5,1.5,4.5,7.5,10.5,13.5]
interval = [-220,-180,-140,-100,-80,-60,-40,-20,0,20,40,60,80,100,140,180,220]
# plot 
diff = np.array(diff)
diff, lon = add_cyclic_point(diff, coord=lon)

# seperate colorbar   
vmax = max(abs(np.reshape(diff,192*289)))
vmin = -max(abs(np.reshape(diff,192*289)))
norm = colors.Normalize(vmin=-250, vmax=250)   

# add figure,coastline&borders
fig = plt.figure(figsize = (24,24))
proj = cartopy.crs.PlateCarree(central_longitude=0)   
ax1 = fig.add_subplot(1, 1, 1, projection=proj)
ax1.outline_patch.set_linewidth(4)
ax1.add_feature(cartopy.feature.COASTLINE,linewidth=4)
reader = Reader('/public/home/gaojy/plot/Data_ipynb/bou1_4p.shp')
provinces = cartopy.feature.ShapelyFeature(
    reader.geometries(), proj, edgecolor='k', facecolor='none')
ax1.add_feature(provinces, linewidth=4)
# set extent
extent=[65,145,10,53]
ax1.set_extent(extent)
# set lon&lat
ax1.set_xticks([75,90,105,120,135])
ax1.set_yticks([15,25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.tick_params(labelsize = 40)


# set title
ax1.set_title(' ',fontdict = {'fontsize' : 45})  

# plot contour
h1 = plt.contourf(lon, lat, diff,
                  interval,
                  cmap = cmaps.MPL_RdBu_r,
                  norm = norm,
                  extend = 'both')

# plot scatter
ax1.plot(index_x,index_y, markersize=10,
         marker='.',color='black',linestyle=' ',transform=proj)

# =============================================================================
# # plot wind change
# skip_xy = (slice(None,None,2))
# skip=(slice(None,None,2),slice(None,None,2))
# Q = ax1.quiver(lon[skip_xy],lat[skip_xy], diff_u[skip], diff_v[skip], scale=22)
# ax1.quiverkey(Q, 0.9, 1.02, 0.5, r'$0.5 m/s$', labelpos='E',fontproperties = {'size':30})
# 
# =============================================================================

#for i in range(len(index_uv_x)):
Q = ax1.quiver(lon_uv[::3],lat_uv[::3],diff_u[::3,::3],diff_v[::3,::3],scale=13,width=0.003,minshaft = 1, minlength=0,headwidth = 2)
ax1.quiverkey(Q, 0.9, 1.02, 1, r'$1 m/s$', labelpos='E',fontproperties = {'size':30})
# add colorbar    
cbar1 = fig.add_axes([0.15,0.18,0.72,0.02]) 
cb1 = plt.colorbar(h1, 
                   cax=cbar1,
                   ticks=interval[::2],
                   orientation='horizontal')
cb1.ax.tick_params(labelsize=40)  
cb1.outline.set_linewidth(2)
cb1.set_label('Pa',fontsize=40)

# save pic
plt.savefig('/public/home/gaojy/plot/CN-PSL-850wind.jpeg',
            #dpi=300, 
            bbox_inches = 'tight')

#plt.show()
