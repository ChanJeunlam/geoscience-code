import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as mticker
import cartopy.feature
from cartopy.util import add_cyclic_point
from matplotlib import colors

emiss1 = np.zeros([2040,13,192,288])
nc1 = nc.Dataset('/public/home/gaojy/CESM2_MEIC_1850_2017/so4_a1_anthro-ene_vertical_mol.nc')
altitude_int = nc1.variables['altitude_int'][:]               
lon1 = nc1.variables['lon'][:]
lat1 = nc1.variables['lat'][:]
for i in range(len(list(nc1.variables.keys()))):
    if (nc1.variables[list(nc1.variables.keys())[i]]).ndim == 4:
        emiss1 += nc1.variables[list(nc1.variables.keys())[i]][:]
emiss1_ = np.mean(emiss1[:,:,:,:] ,0)
emiss1__ = np.zeros([192,288])
for j in range(13):
    emiss1__ += (altitude_int[j+1]-altitude_int[j])*emiss1_[j]*1000000000
emiss1__, lon1 = add_cyclic_point(emiss1__, coord=lon1)

emiss2 = np.zeros([2040,13,192,288])
nc2 = nc.Dataset('/public/home/gaojy/CESM2_MEIC_1850_2017/mod_so4_a1_anthro-ene_vertical_mol.nc')
altitude_int = nc2.variables['altitude_int'][:]               
lon2 = nc2.variables['lon'][:]
lat2 = nc2.variables['lat'][:]
for i in range(len(list(nc2.variables.keys()))):
    if (nc2.variables[list(nc2.variables.keys())[i]]).ndim == 4:
        emiss2 += nc2.variables[list(nc2.variables.keys())[i]][:]
emiss2_ = np.mean(emiss2[:,:,:,:] ,0)
emiss2__ = np.zeros([192,288])
for j in range(13):
    emiss2__ += (altitude_int[j+1]-altitude_int[j])*emiss2_[j]*1000000000
emiss2__, lon2 = add_cyclic_point(emiss2__, coord=lon2)

vmin1 = -max(max(np.reshape(emiss1__,(192*289))), 
           max(np.reshape(emiss2__,(192*289))))
vmax1 = max(max(np.reshape(emiss1__,(192*289))), 
           max(np.reshape(emiss2__,(192*289))))
vmin2 = -max(abs(np.reshape((emiss2__-emiss1__),(192*289))))
vmax2 = max(abs(np.reshape((emiss2__-emiss1__),(192*289))))
norm1 = colors.Normalize(vmin=vmin1, vmax=vmax1) 
norm2 = colors.Normalize(vmin=vmin2, vmax=vmax2) 
  
proj = ccrs.PlateCarree(central_longitude=0)   
fig = plt.figure(figsize = (18,27))
ax1 = fig.add_subplot(3, 1, 1, projection=proj)
ax1.add_feature(cartopy.feature.COASTLINE)
ax1.add_feature(cartopy.feature.BORDERS)
ax1.set_xticks([-160, -120, -80, -40, 0, 40, 80,120, 160])
ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.set(title='so4_a1_anthro-ene_vertical_mol.nc')
h1 = plt.pcolormesh(lon1, lat1, emiss1__, cmap = plt.cm.bwr, norm = norm1)

ax2 = fig.add_subplot(3, 1, 2, projection=proj)
ax2.add_feature(cartopy.feature.COASTLINE)
ax2.add_feature(cartopy.feature.BORDERS)
ax2.set_xticks([-160, -120, -80, -40, 0, 40, 80,120, 160])
ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.set(title='mod_so4_a1_anthro-ene_vertical_mol.nc')
h2 = plt.pcolormesh(lon2, lat2, emiss2__, cmap = plt.cm.bwr, norm = norm1)

ax3 = fig.add_subplot(3, 1, 3, projection=proj)
ax3.add_feature(cartopy.feature.COASTLINE)
ax3.add_feature(cartopy.feature.BORDERS)    
ax3.set_xticks([-160, -120, -80, -40, 0, 40, 80,120, 160])
ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.set(title='latter-former')
h3 = plt.pcolormesh(lon2, lat2, emiss2__-emiss1__, cmap = plt.cm.bwr, 
                    norm = norm2)

#fig.subplots_adjust(right=0.9)

font = {'family' : 'serif',
#       'color'  : 'darkred',
    'color'  : 'black',
    'weight' : 'normal',
    'size'   : 16,
    }
#对应 l,b,w,h；设置colorbar位置；
cbar1 = fig.add_axes([0.92,0.40,0.02,0.48]) 
cb1 = plt.colorbar(h1, cax=cbar1)
cb1.ax.tick_params(labelsize=16)  
cb1.set_label('ppb',fontdict=font) 

cbar2 = fig.add_axes([0.92,0.12,0.02,0.23]) 
cb2 = plt.colorbar(h3, cax=cbar2)
cb2.ax.tick_params(labelsize=16)  
cb2.set_label('ppb',fontdict=font) 

plt.savefig('/public/home/gaojy/CESM2_MEIC_1850_2017/so4_a1_anthro-ene_vertical_mol.png',  dpi=300)
plt.show()

