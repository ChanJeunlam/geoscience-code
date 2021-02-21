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

# define var to be draw 


# put all data into list
base_MAM = []
base_JJA = []
base_SON = []
base_DJF = []
aclean_MAM = []
aclean_JJA = []
aclean_SON = []
aclean_DJF = []
yr = [
      '0006','0007','0008','0009','0010',
      '0011','0012','0013','0014','0015',
      '0016','0017','0018','0019','0020'
     ]
MAM = ['03','04','05']
JJA = ['06','07','08']
SON = ['09','10','11']
DJF = ['01','02','12']
for a in range(len(yr)):
    for b in range(3):
        base_MAM.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[a]+'-'+MAM[b]+'.nc'))
        aclean_MAM.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[a]+'-'+MAM[b]+'.nc'))
        base_JJA.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[a]+'-'+JJA[b]+'.nc'))
        aclean_JJA.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[a]+'-'+JJA[b]+'.nc'))
        base_SON.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[a]+'-'+SON[b]+'.nc'))
        aclean_SON.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[a]+'-'+SON[b]+'.nc'))
        base_DJF.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/Base/atm/hist/Base.cam.h0.'+yr[a]+'-'+DJF[b]+'.nc'))
        aclean_DJF.append(xr.open_dataset
                     ('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.'+yr[a]+'-'+DJF[b]+'.nc'))


aclean_MAM_tdata = []
base_MAM_tdata = []
diff_MAM = np.zeros((192,288))

aclean_JJA_tdata = []
base_JJA_tdata = []
diff_JJA = np.zeros((192,288))

aclean_SON_tdata = []
base_SON_tdata = []
diff_SON = np.zeros((192,288))

aclean_DJF_tdata = []
base_DJF_tdata = []
diff_DJF = np.zeros((192,288)) 

for c in range(len(yr)*3):
    aclean_MAM_ = pd.DataFrame(np.reshape(aclean_MAM[c]['FSNT'].data,(192,288))-np.reshape(aclean_MAM[c]['FSNTC'].data,(192,288))-np.reshape(aclean_MAM[c]['FLNT'].data,(192,288))+np.reshape(aclean_MAM[c]['FLNTC'].data,(192,288))).fillna(0)
    base_MAM_ = pd.DataFrame(np.reshape(base_MAM[c]['FSNT'].data,(192,288))-np.reshape(base_MAM[c]['FSNTC'].data,(192,288))-np.reshape(base_MAM[c]['FLNT'].data,(192,288))+np.reshape(base_MAM[c]['FLNTC'].data,(192,288))).fillna(0) 
    aclean_MAM_tdata.append(aclean_MAM_)
    base_MAM_tdata.append(base_MAM_)
    diff_MAM += (aclean_MAM_-base_MAM_)/(len(yr)*3)
    
    aclean_JJA_ = pd.DataFrame(np.reshape(aclean_JJA[c]['FSNT'].data,(192,288))-np.reshape(aclean_JJA[c]['FSNTC'].data,(192,288))-np.reshape(aclean_JJA[c]['FLNT'].data,(192,288))+np.reshape(aclean_JJA[c]['FLNTC'].data,(192,288))).fillna(0)
    base_JJA_ = pd.DataFrame(np.reshape(base_JJA[c]['FSNT'].data,(192,288))-np.reshape(base_JJA[c]['FSNTC'].data,(192,288))-np.reshape(base_JJA[c]['FLNT'].data,(192,288))+np.reshape(base_JJA[c]['FLNTC'].data,(192,288))).fillna(0) 
    aclean_JJA_tdata.append(aclean_JJA_)
    base_JJA_tdata.append(base_JJA_)
    diff_JJA += (aclean_JJA_-base_JJA_)/(len(yr)*3)

    aclean_SON_ = pd.DataFrame(np.reshape(aclean_SON[c]['FSNT'].data,(192,288))-np.reshape(aclean_SON[c]['FSNTC'].data,(192,288))-np.reshape(aclean_SON[c]['FLNT'].data,(192,288))+np.reshape(aclean_SON[c]['FLNTC'].data,(192,288))).fillna(0)
    base_SON_ = pd.DataFrame(np.reshape(base_SON[c]['FSNT'].data,(192,288))-np.reshape(base_SON[c]['FSNTC'].data,(192,288))-np.reshape(base_SON[c]['FLNT'].data,(192,288))+np.reshape(base_SON[c]['FLNTC'].data,(192,288))).fillna(0) 
    aclean_SON_tdata.append(aclean_SON_)
    base_SON_tdata.append(base_SON_)
    diff_SON += (aclean_SON_-base_SON_)/(len(yr)*3)

    aclean_DJF_ = pd.DataFrame(np.reshape(aclean_DJF[c]['FSNT'].data,(192,288))-np.reshape(aclean_DJF[c]['FSNTC'].data,(192,288))-np.reshape(aclean_DJF[c]['FLNT'].data,(192,288))+np.reshape(aclean_DJF[c]['FLNTC'].data,(192,288))).fillna(0)
    base_DJF_ = pd.DataFrame(np.reshape(base_DJF[c]['FSNT'].data,(192,288))-np.reshape(base_DJF[c]['FSNTC'].data,(192,288))-np.reshape(base_DJF[c]['FLNT'].data,(192,288))+np.reshape(base_DJF[c]['FLNTC'].data,(192,288))).fillna(0) 
    aclean_DJF_tdata.append(aclean_DJF_)
    base_DJF_tdata.append(base_DJF_)
    diff_DJF += (aclean_DJF_-base_DJF_)/(len(yr)*3)
    

aclean_MAM_tdata_mean = ((np.reshape(np.array(aclean_MAM_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)
base_MAM_tdata_mean = ((np.reshape(np.array(base_MAM_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)  

aclean_JJA_tdata_mean = ((np.reshape(np.array(aclean_JJA_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)
base_JJA_tdata_mean = ((np.reshape(np.array(base_JJA_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)

aclean_SON_tdata_mean = ((np.reshape(np.array(aclean_SON_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)
base_SON_tdata_mean = ((np.reshape(np.array(base_SON_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)

aclean_DJF_tdata_mean = ((np.reshape(np.array(aclean_DJF_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)
base_DJF_tdata_mean = ((np.reshape(np.array(base_DJF_tdata),(15,3,192,288))).swapaxes(0,1)).mean(axis=0)




 

lon = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
stat_MAM = np.zeros((192,288))
p_value_MAM = np.zeros((192,288))
index_scatter_MAM = [[],[]]
for p in range(192):
    for q in range(288):    
            stat_MAM[p,q],p_value_MAM[p,q] = stats.ttest_ind(
            base_MAM_tdata_mean[:,p,q],
            aclean_MAM_tdata_mean[:,p,q],
            equal_var=False                   
            )
            if p_value_MAM[p,q]<=0.1:
                index_scatter_MAM[1].append(lat[p])
                index_scatter_MAM[0].append(lon[q])                  
index_scatter_MAM[0].append(-180)
index_scatter_MAM[1].append(90)
index_scatter_MAM[0].append(-180)
index_scatter_MAM[1].append(-90)
index_scatter_MAM[0].append(180)
index_scatter_MAM[1].append(90)
index_scatter_MAM[0].append(180)
index_scatter_MAM[1].append(-90)

stat_JJA = np.zeros((192,288))
p_value_JJA = np.zeros((192,288))
index_scatter_JJA = [[],[]]
for p in range(192):
    for q in range(288):    
            stat_JJA[p,q],p_value_JJA[p,q] = stats.ttest_ind(
            base_JJA_tdata_mean[:,p,q],
            aclean_JJA_tdata_mean[:,p,q],
            equal_var=False                   
            )
            if p_value_JJA[p,q]<=0.1:
                index_scatter_JJA[1].append(lat[p])
                index_scatter_JJA[0].append(lon[q])                  
index_scatter_JJA[0].append(-180)
index_scatter_JJA[1].append(90)
index_scatter_JJA[0].append(-180)
index_scatter_JJA[1].append(-90)
index_scatter_JJA[0].append(180)
index_scatter_JJA[1].append(90)
index_scatter_JJA[0].append(180)
index_scatter_JJA[1].append(-90)

stat_SON = np.zeros((192,288))
p_value_SON = np.zeros((192,288))
index_scatter_SON = [[],[]]
for p in range(192):
    for q in range(288):    
            stat_SON[p,q],p_value_SON[p,q] = stats.ttest_ind(
            base_SON_tdata_mean[:,p,q],
            aclean_SON_tdata_mean[:,p,q],
            equal_var=False                   
            )
            if p_value_SON[p,q]<=0.1:
                index_scatter_SON[1].append(lat[p])
                index_scatter_SON[0].append(lon[q])                  
index_scatter_SON[0].append(-180)
index_scatter_SON[1].append(90)
index_scatter_SON[0].append(-180)
index_scatter_SON[1].append(-90)
index_scatter_SON[0].append(180)
index_scatter_SON[1].append(90)
index_scatter_SON[0].append(180)
index_scatter_SON[1].append(-90)

stat_DJF = np.zeros((192,288))
p_value_DJF = np.zeros((192,288))
index_scatter_DJF = [[],[]]
for p in range(192):
    for q in range(288):    
            stat_DJF[p,q],p_value_DJF[p,q] = stats.ttest_ind(
            base_DJF_tdata_mean[:,p,q],
            aclean_DJF_tdata_mean[:,p,q],
            equal_var=False                   
            )
            if p_value_DJF[p,q]<=0.1:
                index_scatter_DJF[1].append(lat[p])
                index_scatter_DJF[0].append(lon[q])                  
index_scatter_DJF[0].append(-180)
index_scatter_DJF[1].append(90)
index_scatter_DJF[0].append(-180)
index_scatter_DJF[1].append(-90)
index_scatter_DJF[0].append(180)
index_scatter_DJF[1].append(90)
index_scatter_DJF[0].append(180)
index_scatter_DJF[1].append(-90)



lon_MAM = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lon_JJA = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lon_SON = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lon_DJF = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lon.data
lat = (xr.open_dataset('/public/home/gaojy/output/archive/AClean/atm/hist/AClean.cam.h0.0001-01.nc')).lat.data
diff_MAM = np.array(diff_MAM)
diff_MAM, lon_MAM = add_cyclic_point(diff_MAM, coord=lon)
diff_JJA = np.array(diff_JJA)
diff_JJA, lon_JJA = add_cyclic_point(diff_JJA, coord=lon)
diff_SON = np.array(diff_SON)
diff_SON, lon_SON = add_cyclic_point(diff_SON, coord=lon)
diff_DJF = np.array(diff_DJF)
diff_DJF, lon_DJF = add_cyclic_point(diff_DJF, coord=lon)


   
vmax = max(max(abs(np.reshape(diff_MAM,192*289))),max(abs(np.reshape(diff_JJA,192*289)))
           ,max(abs(np.reshape(diff_SON,192*289))),max(abs(np.reshape(diff_DJF,192*289))))
print(vmax)
vmin = -vmax
norm = colors.Normalize(vmin=-9.5, vmax=9.5)   

fig = plt.figure(figsize = (48,30))
proj = cartopy.crs.PlateCarree(central_longitude=0) 
extent=[65,145,10,53]
reader = Reader('/public/home/gaojy/plot/Data_ipynb/bou2_4p.shp')
provinces = cartopy.feature.ShapelyFeature(
    reader.geometries(), proj, edgecolor='k', facecolor='none')
interval = [-6.5,-5.5,-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5]
ax1 = fig.add_subplot(2, 2, 1, projection=proj)
ax1.outline_patch.set_linewidth(2.5)
ax1.add_feature(cartopy.feature.COASTLINE,linewidth=2)
ax1.add_feature(provinces, linewidth=2)
ax1.set_extent(extent)
ax1.set_xticks([75,90,105,120,135])
ax1.set_yticks([15,25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.tick_params(labelsize = 30)
ax1.set_title('MAM',fontdict = {'fontsize' : 30})  
h1 = plt.contourf(lon_MAM, lat, diff_MAM,
                  interval,
                  cmap = plt.cm.bwr,
                  norm = norm,
                  extend = 'both')
ax1.plot(index_scatter_MAM[0],index_scatter_MAM[1], markersize=5,
         marker='.',color='black',linestyle=' ',transform=proj)

ax2 = fig.add_subplot(2, 2, 2, projection=proj)
ax2.outline_patch.set_linewidth(2.5)
ax2.add_feature(cartopy.feature.COASTLINE,linewidth=2)
ax2.add_feature(provinces, linewidth=2)
ax2.set_extent(extent)
ax2.set_xticks([75,90,105,120,135])
ax2.set_yticks([15,25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.tick_params(labelsize =30)
ax2.set_title('JJA',fontdict = {'fontsize' : 30})  
h2 = plt.contourf(lon_JJA, lat, diff_JJA,
                  interval,
                  cmap = plt.cm.bwr,
                  norm = norm,
                  extend = 'both')
ax2.plot(index_scatter_JJA[0],index_scatter_JJA[1], markersize=5,
         marker='.',color='black',linestyle=' ',transform=proj)

ax3 = fig.add_subplot(2, 2, 3, projection=proj)
ax3.outline_patch.set_linewidth(2.5)
ax3.add_feature(cartopy.feature.COASTLINE,linewidth=2)
ax3.add_feature(provinces, linewidth=(2))
ax3.set_extent(extent)
ax3.set_xticks([75,90,105,120,135])
ax3.set_yticks([15,25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.tick_params(labelsize = 30)
ax3.set_title('SON',fontdict = {'fontsize' : 30})  
h3 = plt.contourf(lon_SON, lat, diff_SON,
                  interval,
                  cmap = plt.cm.bwr,
                  norm = norm,
                  extend = 'both')
ax3.plot(index_scatter_SON[0],index_scatter_SON[1], markersize=5,
         marker='.',color='black',linestyle=' ',transform=proj)

ax4 = fig.add_subplot(2, 2, 4, projection=proj)
ax4.outline_patch.set_linewidth(2.5)
ax4.add_feature(cartopy.feature.COASTLINE,linewidth=2)
ax4.add_feature(provinces, linewidth=2)
ax4.set_extent(extent)
ax4.set_xticks([75,90,105,120,135])
ax4.set_yticks([15,25,35,45,55])
lon_formatter = mticker.LongitudeFormatter()
lat_formatter = mticker.LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.tick_params(labelsize = 30)
ax4.set_title('DJF',fontdict = {'fontsize' : 30})  
h4 = plt.contourf(lon_DJF, lat, diff_DJF,
                  interval,
                  cmap = plt.cm.bwr,
                  norm = norm,
                  extend = 'both')
ax4.plot(index_scatter_DJF[0],index_scatter_DJF[1], markersize=5,
         marker='.',color='black',linestyle=' ',transform=proj)




cbar1 = fig.add_axes([0.15,0.00,0.72,0.02]) 
cb1 = plt.colorbar(h1, 
                   cax=cbar1,
                   ticks=interval,
                   orientation='horizontal')
cb1.ax.tick_params(labelsize=30)  
cb1.outline.set_linewidth(2)
cb1.set_label('w m⁻²',fontsize=28)

plt.savefig('/public/home/gaojy/plot/CN-diff-FSNTC-FSNTC_d1-season.jpeg',
            #dpi=300, 
            bbox_inches = 'tight'
            )

#plt.show()
