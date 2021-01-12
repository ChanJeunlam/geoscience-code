import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as mticker
import cartopy.feature
from cartopy.util import add_cyclic_point
from matplotlib import colors
import os
filename3 = []
filename1 = os.listdir("F:/MEIC")
filename2 = [n for n in filename1  if 'anthro' in n and 'vertical' not in n]
for j in range(len(filename2)):
    filename3.append(os.path.splitext(filename2[j])[0])
for k in range(len(filename3)):
    emiss1 = np.zeros([2040,192,288])
    nc1 = nc.Dataset('F:/MEIC/'+filename3[k]+'.nc')               
    lon1 = nc1.variables['lon'][:]
    lat1 = nc1.variables['lat'][:]
    for i in range(len(list(nc1.variables.keys()))):
        if (nc1.variables[list(nc1.variables.keys())[i]]).ndim == 3:
            emiss1 += nc1.variables[list(nc1.variables.keys())[i]][:]
    emiss1_ = np.mean(emiss1[:,:,:] ,0)
    emiss1_, lon1 = add_cyclic_point(emiss1_, coord=lon1)
    
    emiss2 = np.zeros([2040,192,288])
    nc2 = nc.Dataset('F:/mod/mod_'+filename3[k]+'.nc')
    lon2 = nc2.variables['lon'][:]
    lat2 = nc2.variables['lat'][:]
    for i in range(len(list(nc2.variables.keys()))):
        if (nc2.variables[list(nc2.variables.keys())[i]]).ndim == 3:
            emiss2 += nc2.variables[list(nc2.variables.keys())[i]][:]
    emiss2_ = np.mean(emiss2[:,:,:] ,0)  
    emiss2_, lon2 = add_cyclic_point(emiss2_, coord=lon2)
    
    vmin1 = -max(max(np.reshape(emiss1_,(192*289))), 
               max(np.reshape(emiss2_,(192*289))))
    vmax1 = max(max(np.reshape(emiss1_,(192*289))), 
               max(np.reshape(emiss2_,(192*289))))
    vmin2 = -max(abs(np.reshape((emiss2_-emiss1_),(192*289))))
    vmax2 = max(abs(np.reshape((emiss2_-emiss1_),(192*289))))
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
    ax1.set(title=filename3[k])
    h1 = plt.pcolormesh(lon1, lat1, emiss1_, cmap = plt.cm.bwr, norm = norm1)
    
    ax2 = fig.add_subplot(3, 1, 2, projection=proj)
    ax2.add_feature(cartopy.feature.COASTLINE)
    ax2.add_feature(cartopy.feature.BORDERS)
    ax2.set_xticks([-160, -120, -80, -40, 0, 40, 80,120, 160])
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    lon_formatter = mticker.LongitudeFormatter()
    lat_formatter = mticker.LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    ax2.set(title='mod_'+filename3[k])
    h2 = plt.pcolormesh(lon2, lat2, emiss2_, cmap = plt.cm.bwr, norm = norm1)
    
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
    h3 = plt.pcolormesh(lon2, lat2, emiss2_-emiss1_, cmap = plt.cm.bwr, 
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
    cb1.set_label('molecules/cm²/s',fontdict=font) 
    
    cbar2 = fig.add_axes([0.92,0.12,0.02,0.23]) 
    cb2 = plt.colorbar(h3, cax=cbar2)
    cb2.ax.tick_params(labelsize=16)  
    cb2.set_label('molecules/cm²/s',fontdict=font) 
    
    plt.savefig('F:/MEIC/'+filename3[k]+'.png',  dpi=300)
    plt.show()
    
