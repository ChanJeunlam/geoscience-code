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

nc2013 = xr.open_dataset('F:/2013_so4_a1_anthro-ag-ship_surface_mol.nc')
nc2017 = xr.open_dataset('F:/2017_mod_so4_a1_anthro-ag-ship_surface_mol.nc')
var2013 = list(nc2013.data_vars)
var2017 = list(nc2017.data_vars)
nc2013_ = xr.open_dataset('F:/2013_so4_a2_anthro-res_surface_mol.nc')
nc2017_ = xr.open_dataset('F:/2017_mod_so4_a2_anthro-res_surface_mol.nc')
var2013_ = list(nc2013_.data_vars)
var2017_= list(nc2017_.data_vars)
nc2013__ = xr.open_dataset('F:/2013_so4_a1_anthro-ene_vertical_mol.nc')
nc2017__ = xr.open_dataset('F:/2017_mod_so4_a1_anthro-ene_vertical_mol.nc')
var2013__ = list(nc2013__.data_vars)
var2017__ = list(nc2017__.data_vars)

data2013 = np.zeros((192,288))
data2017 = np.zeros((192,288))
for i in range(0,len(var2013)): 
    if nc2013[var2013[i]].dims == ('time', 'lat', 'lon'):    
        data2013 += np.mean(nc2013[var2013[i]].data,0)
for i in range(0,len(var2013)): 
    if nc2017[var2017[i]].dims == ('time', 'lat', 'lon'):    
        data2017 += np.mean(nc2017[var2017[i]].data,0)
        
for i in range(0,len(var2013_)): 
    if nc2013_[var2013_[i]].dims == ('time', 'lat', 'lon'):    
        data2013 += np.mean(nc2013_[var2013_[i]].data,0)
for i in range(0,len(var2013_)): 
    if nc2017_[var2017_[i]].dims == ('time', 'lat', 'lon'):    
        data2017 += np.mean(nc2017_[var2017_[i]].data,0)
        
for i in range(0,len(var2013__)): 
    if nc2013__[var2013__[i]].dims == ('time', 'altitude','lat', 'lon'):    
        data2013 += np.mean(nc2013__[var2013__[i]].data,0)[0]
for i in range(0,len(var2013__)): 
    if nc2017__[var2017__[i]].dims == ('time','altitude', 'lat', 'lon'):    
        data2017 += np.mean(nc2017__[var2017__[i]].data,0)[0]

lon = (xr.open_dataset('F:/2013_bc_a4_anthro_surface_mol.nc')).lon.data
lat = (xr.open_dataset('F:/2013_bc_a4_anthro_surface_mol.nc')).lat.data


diff = (data2017-data2013)/6.022e23*115*1e3*1e4*60*60*24

# plot 
diff = np.array(diff)
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
ax1.set_title('(d)SO₄²⁻',fontdict = {'fontsize' : 20},loc = 'left')  
# seperate colorbar   
vmax = max(abs(np.reshape(diff,192*289)))
vmin = -max(abs(np.reshape(diff,192*289)))
norm = colors.Normalize(vmin=-0.5, vmax=0.5) 
#interval = [-1.35,-1.05,-0.75,-0.45,-0.15,-0.09,-0.03,0.03,0.09,0.15,0.45,0.75,1.05,1.35]
#interval = [-2.25,-1.95,-1.65,-1.35,-1.05,-0.75,-0.45,-0.15,0.15]
#interval = [-5.2,-4.4,-3.6,-2.8,-2.0,-1.2,-0.4,0.4]
#interval = [-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1]
interval = [-0.45,-0.39,-0.33,-0.27,-0.21,-0.15,-0.09,-0.03,0.03]
# plot contour
h1 = plt.contourf(lon, lat, diff,
                  interval,
                  cmap = cmaps.MPL_bwr,
                  norm = norm,
                  #extend = 'min'
                  )
#add colorba
cbar1 = fig.add_axes([0.15,0.18,0.72,0.02]) 
cb1 = plt.colorbar(h1, 
                   cax=cbar1,
                   ticks=interval,
                   orientation='horizontal')
cb1.ax.tick_params(labelsize=13)  
cb1.outline.set_linewidth(1)
cb1.set_label('mg m⁻² day⁻¹',fontsize=15)
# save pic
plt.savefig('E:/so4.jpeg',
            #dpi=300, 
            bbox_inches = 'tight')

#plt.show()
