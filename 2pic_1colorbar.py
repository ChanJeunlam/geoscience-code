import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as mticker
import cartopy.feature
from cartopy.util import add_cyclic_point
from matplotlib import colors


emiss1 = np.zeros([12,192,288])
nc1 = nc.Dataset('F:/SOAGx1.5_anthro_surface_mol.nc')               
lon1 = nc1.variables['lon'][:]
lat1 = nc1.variables['lat'][:]
for i in range(len(list(nc1.variables.keys()))):
    if (nc1.variables[list(nc1.variables.keys())[i]]).ndim == 3:
        emiss1 += nc1.variables[list(nc1.variables.keys())[i]][1968:1980]
emiss1_ = np.mean(emiss1[:,:,:] ,0)
emiss1_, lon1 = add_cyclic_point(emiss1_, coord=lon1)

emiss2 = np.zeros([12,192,288])
nc2 = nc.Dataset('F:/r_mod_SOAGx1.5_anthro_surface_mol.nc')
lon2 = nc2.variables['lon'][:]
lat2 = nc2.variables['lat'][:]
for i in range(len(list(nc2.variables.keys()))):
    if (nc2.variables[list(nc2.variables.keys())[i]]).ndim == 3:
        emiss2 += nc2.variables[list(nc2.variables.keys())[i]][1968:1980]
emiss2_ = np.mean(emiss2[:,:,:] ,0)  
emiss2_, lon2 = add_cyclic_point(emiss2_, coord=lon2)

emiss3 = np.zeros([12,192,288])
nc3 = nc.Dataset('F:/SOAG_EA_L.nc')
lon3 = nc3.variables['lon'][:]
lat3 = nc3.variables['lat'][:]
for i in range(len(list(nc3.variables.keys()))):
    if (nc3.variables[list(nc3.variables.keys())[i]]).ndim == 3:
        emiss3 += nc3.variables[list(nc3.variables.keys())[i]][0:12]
emiss3_ = np.mean(emiss3[:,:,:] ,0)  
emiss3_, lon3 = add_cyclic_point(emiss3_, coord=lon3)


emiss_gao = (emiss2_-emiss1_)/10**11
emiss_liu = (emiss3_-emiss1_)/10**11

vmin = -max(max(np.reshape(abs(emiss_gao),(192*289))), 
           max(np.reshape(abs(emiss_liu),(192*289))))
vmax = max(max(np.reshape(abs(emiss_gao),(192*289))), 
           max(np.reshape(abs(emiss_liu),(192*289))))

norm = colors.Normalize(vmin=vmin, vmax=vmax) 
  
proj = ccrs.PlateCarree(central_longitude=0)   
fig = plt.figure(figsize = (10,12))
ax1 = fig.add_subplot(2, 1, 1, projection=proj)
ax1.add_feature(cartopy.feature.COASTLINE)
ax1.add_feature(cartopy.feature.BORDERS)
extent=[70,150,15,55]
ax1.set_extent(extent)
ax1.set_xticks([70,90,110,130,150])
ax1.set_yticks([25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.set(title='gao')
h1 = plt.contourf(lon1, lat1, emiss_gao,
                  [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2],
                  cmap = plt.cm.bwr,
                  norm = norm
                  )
ax2 = fig.add_subplot(2, 1, 2, projection=proj)
ax2.add_feature(cartopy.feature.COASTLINE)
ax2.add_feature(cartopy.feature.BORDERS)
extent=[70,150,15,55]
ax2.set_extent(extent)
ax2.set_xticks([70,90,110,130,150])
ax2.set_yticks([25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.set(title='liu')
h2 = plt.contourf(lon2, lat2, emiss_liu,
                  [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2],
                  cmap = plt.cm.bwr,
                  norm = norm
                  )


#fig.subplots_adjust(right=0.9)

font = {'family' : 'serif',
#       'color'  : 'darkred',
    'color'  : 'black',
    'weight' : 'normal',
    'size'   : 16,
    }

cbar1 = fig.add_axes([0.15,0.075,0.7,0.015]) 
cb1 = plt.colorbar(h1, cax=cbar1,
  ticks=[-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2],
  orientation='horizontal')
cb1.ax.tick_params(labelsize=16)  
cb1.set_label('10¹¹ molecules/cm²/s',fontdict=font) 

'''
cbar2 = fig.add_axes([0.15,0.2,0.7,0.015]) 
cb2 = plt.colorbar(h2, cax=cbar2,orientation='horizontal')
cb2.ax.tick_params(labelsize=16)  
cb2.set_label('molecules/cm²/s',fontdict=font) 
'''


plt.savefig('F:/test_SOAG.png',  dpi=300)
plt.show()

