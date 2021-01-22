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

# put all data into a list
yr = [
      #'0001','0002','0003','0004','0005',
      '0006','0007','0008','0009','0010',
      '0011','0012','0013','0014'
     ]
mon = [
       '01','02','03','04','05','06','07','08','09','10','11','12'
      ]
base = []
for i in range(len(yr)):
    for j in range(len(mon)):
        base.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[i]+'-'+mon[j]+'.nc'))   
aclean = []
for i in range(len(yr)):
    for j in range(len(mon)):
        aclean.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[i]+'-'+mon[j]+'.nc'))
        
        
        
lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data    
     
# calc conc diff(var indexed by m)&data needed in t-test
var = [
       #'bc_a1_SRF','bc_a4_SRF','pom_a1_SRF','pom_a4_SRF','so4_a1_SRF',
       #'so4_a2_SRF','so4_a3_SRF','soa_a1_SRF','soa_a2_SRF'
       'AEROD_v'
      ]
nmon = len(yr)*len(mon)
    # ! attention:data length decided by len(var)
acleandata = [[]]          
basedata = [[]]
conc_diff = np.zeros((len(var),len(lat),len(lon)))
for m in range(len(var)):
    for n in range(nmon):
        aclean_ = pd.DataFrame(np.reshape(aclean[n][var[m]].data,(192,288))).fillna(0) *np.reshape(aclean[n]['PS'].data,(192,288))/np.reshape(aclean[n]['TS'].data,(192,288))/287.05*1e9 
        base_ = pd.DataFrame(np.reshape(base[n][var[m]].data,(192,288))).fillna(0) *np.reshape(base[n]['PS'].data,(192,288))/np.reshape(base[n]['TS'].data,(192,288))/287.05*1e9
        acleandata[m].append(aclean_)
        basedata[m].append(base_)
        conc_diff[m] += (aclean_-base_)/nmon                  
    #data into array
acleandata = np.array(acleandata)
basedata = np.array(basedata)  




# calc p value in t-test      
stat = np.zeros((len(var),len(lat),len(lon)))
p = np.zeros((len(var),len(lat),len(lon)))
    # ! attention:data length decided by len(var)
index_scatter = [[[],[]]]
for a in range(len(var)):
    for b in range(len(lat)):
        for c in range(len(lon)):    
                stat[a,b,c],p[a,b,c] = stats.ttest_ind(
                basedata[a,:,b,c],
                acleandata[a,:,b,c],
                equal_var=False                       )
                if p[a,b,c]<=0.05:
                    index_scatter[a][1].append(lat[b])
                    index_scatter[a][0].append(lon[c]) 
    # add 4 points to extend fig                    
    index_scatter[a][0].append(-180)
    index_scatter[a][1].append(90)
    index_scatter[a][0].append(-180)
    index_scatter[a][1].append(-90)
    index_scatter[a][0].append(180)
    index_scatter[a][1].append(90)
    index_scatter[a][0].append(180)
    index_scatter[a][1].append(-90)



# plot horizontal distribution     
for d in range(len(var)):
    # remove middle white line
    lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
    lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
    conc_ = conc_diff[d]
    conc_, lon = add_cyclic_point(conc_, coord=lon)
    
    # seperate colorbar into plus&minus  
    vmax = max(abs(np.reshape(conc_,192*289)))
    vmin = -max(abs(np.reshape(conc_,192*289)))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)   
    
    # add figure,coastline&borders
    fig = plt.figure(figsize = (24,18))
    proj = cartopy.crs.PlateCarree(central_longitude=0)   
    ax1 = fig.add_subplot(1, 1, 1, projection=proj)
    ax1.outline_patch.set_linewidth(4)
    #ax1.spines['left'].set_color('red')
    #ax1.spines['left'].set_linewidth(2)
    ax1.add_feature(cartopy.feature.COASTLINE,linewidth=4)
    ax1.add_feature(cartopy.feature.BORDERS,linewidth=4)
    
    # set extent
    extent=[70,140,15,50]
    ax1.set_extent(extent)
    
    # set lon&lat
    #ax1.set_xticks([-160,-120,-80,-40,0,40,80,120,160])
    #ax1.set_yticks([-90,-60,-30,0,30,60,90])
    ax1.set_xticks([75,90,105,120,135])
    ax1.set_yticks([20,30,40,50])
    lon_formatter = mticker.LongitudeFormatter()
    lat_formatter = mticker.LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    ax1.tick_params(labelsize = 40)
    
    # set title
    ax1.set_title(var[d],fontdict = {'fontsize' : 45})  
    
    # set contour interval
    interval = [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2]
    
    # plot contour
    h1 = plt.contourf(lon, lat, conc_,
                      #interval,
                      cmap = plt.cm.bwr,
                      norm = norm)
    # plot scatter
    ax1.plot(index_scatter[d][0],index_scatter[d][1], markersize=10,
             marker='.',color='black',linestyle=' ',transform=proj)
                         
    # add colorbar    
    cbar1 = fig.add_axes([0.15,0.12,0.72,0.03]) 
    cb1 = plt.colorbar(h1, 
                       cax=cbar1,
                       #ticks=interval,
                       orientation='horizontal')
    cb1.ax.tick_params(labelsize=40)  
    cb1.outline.set_linewidth(2)
    cb1.set_label(' ',fontsize=40)
    
    # save fig
    plt.savefig('/public/home/gaojy/plot/CN_conc_diff_'+var[d]+'.png',
                #dpi=300, 
                bbox_inches = 'tight')
    #plt.show()

