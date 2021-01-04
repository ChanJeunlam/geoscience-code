import xarray as xr
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import cartopy
import cartopy.mpl.ticker as mticker
mon = ['01','02','03','04','05','06','07','08','09','10','11','12']

conc_Base = []
for i in range(len(mon)):
    conc_Base.append(xr.open_dataset('F:/Base.cam.h0.0001-'+mon[i]+'.nc'))
    
conc_AClean = []
for i in range(len(mon)):
    conc_AClean.append(xr.open_dataset('F:/AClean.cam.h0.0001-'+mon[i]+'.nc'))

var = ['bc_a1_SRF','bc_a4_SRF','pom_a1_SRF','pom_a4_SRF','so4_a1_SRF',
       'soa_a1_SRF','AEROD_v']

conc_diff = np.zeros((7,192,288))
for j in range(len(var)):
    for i in range(len(mon)):
        conc_diff[j] +=  (np.reshape(conc_AClean[i][var[j]].data,(192,288))
                          -np.reshape(conc_Base[i][var[j]].data,(192,288)))
                    


    
for j in range(len(var)):
    lon = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lon.data
    lat = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lat.data
#    conc_diff[j], lon = add_cyclic_point(conc_diff[j], coord=lon)
    vmax = max(abs(np.reshape(conc_diff[j],192*288)))
    vmin = -max(abs(np.reshape(conc_diff[j],192*288)))

    norm = colors.Normalize(vmin=vmin, vmax=vmax)   
    proj = cartopy.crs.PlateCarree(central_longitude=0)   
    fig = plt.figure(figsize = (10,9))
    ax1 = fig.add_subplot(1, 1, 1, projection=proj)
    ax1.add_feature(cartopy.feature.COASTLINE)
    ax1.add_feature(cartopy.feature.BORDERS)
#    extent=[70,150,15,55]
#    ax1.set_extent(extent)
    ax1.set_xticks([70,90,110,130,150])
    ax1.set_yticks([25,35,45,55])
    lon_formatter = mticker.LongitudeFormatter()
    lat_formatter = mticker.LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    ax1.set(title=var[j])
    h1 = plt.contourf(lon, lat, conc_diff[j],
#                      [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2],
                      cmap = plt.cm.bwr,
                      norm = norm
                      )
                         
    

    font = {'family' : 'serif',
    #       'color'  : 'darkred',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }
    
    cbar1 = fig.add_axes([0.15,0.075,0.7,0.015]) 
    cb1 = plt.colorbar(h1, cax=cbar1,
    #  ticks=[-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2],
      orientation='horizontal')
    cb1.ax.tick_params(labelsize=16)  
    cb1.set_label('  ',fontdict=font)
     
    #plt.savefig('F:/test.png',  dpi=300)
    plt.show()