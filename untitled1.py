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
# add figure,coastline&borders
fig = plt.figure(figsize = (8,6))
proj = cartopy.crs.PlateCarree(central_longitude=0)   
ax1 = fig.add_subplot(1, 1, 1, projection=proj)
x = [-180,-180,180,180,0]
y = [90,-90,90,-90,-60]

ax1.coastlines()
#ax1.add_feature(cartopy.feature.COASTLINE)
#ax1.add_feature(cartopy.feature.BORDERS)
ax1.scatter(x, y,s = 2,color = 'black'
            #markersize=0.1,marker='o',linestyle='',color='black'
            )

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
#ax1.set_title(v,fontdict = {'fontsize' : 18})




plt.show()