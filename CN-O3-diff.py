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


nc2013 = xr.open_dataset('E:/O3_2013_GC4CESM2.nc')
nc2017 = xr.open_dataset('E:/O3_2017_GC4CESM2.nc')

lon = (xr.open_dataset('E:/O3_2017_GC4CESM2.nc')).lon.data
lat = (xr.open_dataset('E:/O3_2017_GC4CESM2.nc')).lat.data
diff = np.mean(nc2017.O3.data-nc2013.O3.data,0)[0]*1e9

diff, lon = add_cyclic_point(diff, coord=lon)
# add figure,coastline&borders
fig = plt.figure(figsize = (12,12))
proj = cartopy.crs.PlateCarree(central_longitude=0)   
ax1 = fig.add_subplot(1, 1, 1, projection=proj)
ax1.outline_patch.set_linewidth(2)
ax1.add_feature(cartopy.feature.COASTLINE,linewidth=1.5)
reader = Reader('E:/Data_ipynb/bou1_4p.shp')
provinces = cartopy.feature.ShapelyFeature(
    reader.geometries(), proj, edgecolor='k', facecolor='none')
ax1.add_feature(provinces, linewidth=1.5)
#ax1.add_feature(cartopy.feature.BORDERS,linewidth=4)

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
ax1.tick_params(labelsize = 20)
# set title
ax1.set_title(' ',fontdict = {'fontsize' : 20},loc = 'left')  
# seperate colorbar   
vmax = max(abs(np.reshape(diff,192*289)))
vmin = -max(abs(np.reshape(diff,192*289)))
norm = colors.Normalize(vmin=vmin, vmax=vmax) 
interval = [-1.35,-1.05,-0.75,-0.45,-0.15,-0.09,-0.03,0.03,0.09,0.15,0.45,0.75,1.05,1.35]
# plot contour
h1 = plt.contourf(lon, lat, diff,
                  #interval,
                  cmap = cmaps.MPL_bwr,
                  norm = norm,
                  #extend = 'both'
                  )
#add colorba
cbar1 = fig.add_axes([0.15,0.18,0.72,0.02]) 
cb1 = plt.colorbar(h1, 
                   cax=cbar1,
                   #ticks=interval,
                   orientation='horizontal')
cb1.ax.tick_params(labelsize=20)  
cb1.outline.set_linewidth(1)
cb1.set_label('ppbv',fontsize=15)
# save pic
plt.savefig('E:/O3_diff.jpeg',
            #dpi=300, 
            bbox_inches = 'tight')

#plt.show()

