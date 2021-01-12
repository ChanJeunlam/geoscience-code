import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netCDF4 as nc

modelID = ['CanESM5','CESM2-WACCM','GFDL-ESM4','INM-CM4-8','INM-CM5-0','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']
sspID = ['ssp126','ssp245','ssp370','ssp585']
ssp_ID = ['ssp126','ssp245','ssp370','ssp585', 'ssp126','ssp245','ssp370','ssp585']
legend = [0, 0, 0, 0, 0, 0, 0, 0]

time = list(range(2015, 2101))
fig, ax = plt.subplots(2, 4,  figsize = (40,20))

for i in range(2):
    for j in range(4):
        nc_126 = nc.Dataset('regrid_emidust_AERmon_'+modelID[4*i+j]+'_'+sspID[0]+'_r1i1p1f1_201501-210012.nc')
        nc_245 = nc.Dataset('regrid_emidust_AERmon_'+modelID[4*i+j]+'_'+sspID[1]+'_r1i1p1f1_201501-210012.nc')
        nc_370 = nc.Dataset('regrid_emidust_AERmon_'+modelID[4*i+j]+'_'+sspID[2]+'_r1i1p1f1_201501-210012.nc')
        nc_585 = nc.Dataset('regrid_emidust_AERmon_'+modelID[4*i+j]+'_'+sspID[3]+'_r1i1p1f1_201501-210012.nc')
        nc_area = nc.Dataset('area_f09.nc')
        emidust126 = nc_126.variables['emidust'][:]*nc_area.variables['area'][:]
        emidust245 = nc_245.variables['emidust'][:]*nc_area.variables['area'][:]
        emidust370 = nc_370.variables['emidust'][:]*nc_area.variables['area'][:]
        emidust585 = nc_585.variables['emidust'][:]*nc_area.variables['area'][:]
        
        emidust_126 = (np.sum(np.sum(np.sum(emidust126[:, 42:96, 0:49],1),1).reshape(86,12),1)+np.sum(np.sum(np.sum(emidust126[:, 42:96, 273:288],1),1).reshape(86,12),1))*30*24*60*60*0.000000000001
        emidust_245 = (np.sum(np.sum(np.sum(emidust245[:, 42:96, 0:49],1),1).reshape(86,12),1)+np.sum(np.sum(np.sum(emidust245[:, 42:96, 273:288],1),1).reshape(86,12),1))*30*24*60*60*0.000000000001
        emidust_370 = (np.sum(np.sum(np.sum(emidust370[:, 42:96, 0:49],1),1).reshape(86,12),1)+np.sum(np.sum(np.sum(emidust370[:, 42:96, 273:288],1),1).reshape(86,12),1))*30*24*60*60*0.000000000001
        emidust_585 = (np.sum(np.sum(np.sum(emidust585[:, 42:96, 0:49],1),1).reshape(86,12),1)+np.sum(np.sum(np.sum(emidust585[:, 42:96, 273:288],1),1).reshape(86,12),1))*30*24*60*60*0.000000000001
        #emidust_126 = np.sum(np.sum(np.sum(emidust126[:, 42:84, 80:257],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_245 = np.sum(np.sum(np.sum(emidust245[:, 42:84, 80:257],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_370 = np.sum(np.sum(np.sum(emidust370[:, 42:84, 80:257],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_585 = np.sum(np.sum(np.sum(emidust585[:, 42:84, 80:257],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_126 = np.sum(np.sum(np.sum(emidust126[:],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_245 = np.sum(np.sum(np.sum(emidust245[:],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_370 = np.sum(np.sum(np.sum(emidust370[:],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        #emidust_585 = np.sum(np.sum(np.sum(emidust585[:],1),1).reshape(86,12),1)*30*24*60*60*0.000000000001
        
        
        legend[0] = ax[i, j].plot(time, emidust_126, color='red', linestyle='dashed')
        legend[1] = ax[i, j].plot(time, emidust_245, color='blue', linestyle='dashed')
        legend[2] = ax[i, j].plot(time, emidust_370, color='green', linestyle='dashed')
        legend[3] = ax[i, j].plot(time, emidust_585, color='orange', linestyle='dashed')
        
        p_126 = np.poly1d(np.polyfit(time, emidust_126, 1))
        legend[4] = ax[i, j].plot(time, p_126(time), color='red', linestyle='solid', linewidth=6)
        p_245 = np.poly1d(np.polyfit(time, emidust_245, 1))
        legend[5] = ax[i, j].plot(time, p_245(time), color='blue', linestyle='solid', linewidth=6)
        p_370 = np.poly1d(np.polyfit(time, emidust_370, 1))
        legend[6] = ax[i, j].plot(time, p_370(time), color='green', linestyle='solid', linewidth=6)
        p_585 = np.poly1d(np.polyfit(time, emidust_585, 1))
        legend[7] = ax[i, j].plot(time, p_585(time), color='orange', linestyle='solid', linewidth=6)
        
        ax[i, j].set_title(modelID[4*i+j], fontsize = 28)
        ax[i,j].tick_params(labelsize = 20)
        
        
fig.suptitle('Dust Emissions Based on 4 Scenarios from 8 Models(South Africa)',fontsize = 38)

fig.legend(legend, labels=ssp_ID, loc = 'upper right', fontsize = 30)

fig.text(0.5, 0.08, 'Year',  fontsize = 28)
fig.text(0.09, 0.4, 'Dust Emissions(Tg/yr)', rotation='vertical', fontsize = 28)

plt.savefig("Dust Emissions Based on 4 Scenarios from 8 Models(South Africa).png")
plt.show()
