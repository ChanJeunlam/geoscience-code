import xarray as xr
import numpy as np
nc1 = xr.open_dataset('F:/pom_a4_anthro_surface_mol.nc')
nc2 = xr.open_dataset('E:/Archive/area_f09.nc')
ENE2013 = []
IND2013 = []
TRA2013 = []
RCO2013 = []
WST2013 = []
SHP2013 = []
for i in range(12):
    ENE2013.append(nc1['ENE'].data[1968+i]*nc2['area'].data)
    IND2013.append(nc1['IND'].data[1968+i]*nc2['area'].data)
    TRA2013.append(nc1['TRA'].data[1968+i]*nc2['area'].data)
    RCO2013.append(nc1['RCO'].data[1968+i]*nc2['area'].data)
    WST2013.append(nc1['WST'].data[1968+i]*nc2['area'].data)
    SHP2013.append(nc1['SHP'].data[1968+i]*nc2['area'].data)
ENE2013_sum = np.sum(np.sum(ENE2013,1),1)
IND2013_sum = np.sum(np.sum(IND2013,1),1)
TRA2013_sum = np.sum(np.sum(TRA2013,1),1)
RCO2013_sum = np.sum(np.sum(RCO2013,1),1)
WST2013_sum = np.sum(np.sum(WST2013,1),1)
SHP2013_sum = np.sum(np.sum(SHP2013,1),1)

ENE2017 = []
IND2017 = []
TRA2017 = []
RCO2017 = []
WST2017 = []
SHP2017 = []
for i in range(12):
    ENE2017.append(nc1['ENE'].data[2016+i]*nc2['area'].data)
    IND2017.append(nc1['IND'].data[2016+i]*nc2['area'].data)
    TRA2017.append(nc1['TRA'].data[2016+i]*nc2['area'].data)
    RCO2017.append(nc1['RCO'].data[2016+i]*nc2['area'].data)
    WST2017.append(nc1['WST'].data[2016+i]*nc2['area'].data)
    SHP2017.append(nc1['SHP'].data[2016+i]*nc2['area'].data)
ENE2017_sum = np.sum(np.sum(ENE2017,1),1)
IND2017_sum = np.sum(np.sum(IND2017,1),1)
TRA2017_sum = np.sum(np.sum(TRA2017,1),1)
RCO2017_sum = np.sum(np.sum(RCO2017,1),1)
WST2017_sum = np.sum(np.sum(WST2017,1),1)
SHP2017_sum = np.sum(np.sum(SHP2017,1),1)

r_ENE = ENE2017_sum/ENE2013_sum
r_IND = IND2017_sum/IND2013_sum
r_TRA = TRA2017_sum/TRA2013_sum
r_RCO = RCO2017_sum/RCO2013_sum
r_WST = WST2017_sum/WST2013_sum
r_SHP = SHP2017_sum/SHP2013_sum

print(r_ENE)
print(r_IND)
print(r_TRA)
print(r_RCO)
print(r_WST)
print(r_SHP)

flag_nc = xr.open_dataset('F:/number_29_tag_regions.nc')
nc4 = xr.open_dataset('F:/SOAGx1.5_anthro_surface_EAN.nc')
nc3 = xr.open_dataset('F:/SOAGx1.5_anthro_surface_mol.nc')
for i in range(192):
    for j in range(288):
        if (flag_nc['region_code'].data[i,j]>=61 
                           and flag_nc['region_code'].data[i,j]<=68):
            for k in range(12):
                nc4['ENE'].data[k,i,j] = nc3['ENE'].data[1968+k,i,j]*r_ENE[k]
                nc4['IND'].data[k,i,j] = nc3['IND'].data[1968+k,i,j]*r_IND[k]
                nc4['TRA'].data[k,i,j] = nc3['TRA'].data[1968+k,i,j]*r_TRA[k]
                nc4['RCO'].data[k,i,j] = nc3['RCO'].data[1968+k,i,j]*r_RCO[k]
                nc4['WST'].data[k,i,j] = nc3['WST'].data[1968+k,i,j]*r_WST[k]
                nc4['SHP'].data[k,i,j] = nc3['SHP'].data[1968+k,i,j]*r_SHP[k]

nc4.to_netcdf('F:/liu_r_mod_SOAGx1.5_anthro_surface_mol.nc','w','NETCDF4_CLASSIC')


