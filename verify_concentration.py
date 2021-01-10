# import packages
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



# calculate concentration differece
    #calculate concentration difference individually
var = [
       'bc_a1_SRF','bc_a4_SRF','pom_a1_SRF','pom_a4_SRF','so4_a1_SRF',
       'so4_a2_SRF','so4_a3_SRF','soa_a1_SRF','soa_a2_SRF'
       #AEROD_v
       ]

conc_diff = np.zeros((10,192,288))
for m in range(len(var)):
    for n in range(len(Base)):
        conc_diff[m] +=  (
            ( pd.DataFrame(np.reshape(AClean[n][var[m]].data,(192,288))).fillna(0) ) *np.reshape(AClean[n]['PS'].data,(192,288))/np.reshape(AClean[n]['TS'].data,(192,288))/287.05*1e9 
           -( pd.DataFrame(np.reshape(Base[n][var[m]].data,(192,288))).fillna(0)) *np.reshape(Base[n]['PS'].data,(192,288))/np.reshape(Base[n]['TS'].data,(192,288))/287.05*1e9
                         )
# =============================================================================
#     #calculate concentration difference totally
# var_ = ['bc_SRF','pom_SRF','so4_SRF','soa_SRF']
# conc_diff_ =[]
# conc_diff_.append(conc_diff[0]+conc_diff[1])
# conc_diff_.append(conc_diff[2]+conc_diff[3])
# conc_diff_.append(conc_diff[4]+conc_diff[5]+conc_diff[6])
# conc_diff_.append(conc_diff[7]+conc_diff[8])   
#                 
# =============================================================================
# plot horizontal distribution     
for m in range(len(var)):
    # remove middle white line
    lon = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lon.data
    lat = (xr.open_dataset('F:/AClean.cam.h0.0001-01.nc')).lat.data
    conc_ = conc_diff[m]
    conc_, lon = add_cyclic_point(conc_, coord=lon)
    
    # make colorbar white in 0
    vmax = max(abs(np.reshape(conc_,192*289)))
    vmin = -max(abs(np.reshape(conc_,192*289)))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)   
    
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
    
    # set contour interval
    interval = [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2]
    
    # plot contour
    h1 = plt.contourf(lon, lat, conc_,
                      #interval,
                      cmap = plt.cm.bwr,
                      norm = norm)
                         
    
    # add colorbar    
    cbar1 = fig.add_axes([0.15,0.12,0.7,0.02]) 
    cb1 = plt.colorbar(h1, 
                       cax=cbar1,
                       #ticks=interval,
                       orientation='horizontal')
    cb1.ax.tick_params(labelsize=10)  
    cb1.set_label('µg/m³')
     
    plt.savefig('F:/conc_diff_'+var[m]+'.png',
                dpi=300, 
                bbox_inches = 'tight')
    plt.show()

