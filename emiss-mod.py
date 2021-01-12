import xarray as xr



flag_nc = xr.open_dataset('E:/number_29_tag_regions.nc')
emiss1 = xr.open_dataset('E:/2013_num_so4_a2_anthro-res_surface_mol.nc')
emiss2 = xr.open_dataset('E:/2017_num_so4_a2_anthro-res_surface_mol.nc')
lat = emiss1['lat'].data
lon = emiss1['lon'].data
sec = ['ENE','IND','TRA','RCO','SLV','WST','SHP']
for k in range(0,len(sec)):
    for i in range(0,len(lat)):
        for j in range(0,len(lon)):
            if (flag_nc['region_code'].data[i,j]>=61 
                and flag_nc['region_code'].data[i,j]<=68):
                if sec[k] in list(emiss1.data_vars):
                    emiss1[sec[k]].data[0:12,i,j] = (emiss2[sec[k]]
                                                     .data[0:12,i,j])
    emiss1.to_netcdf('E:/'+sec[k]+
                     '_AClean_num_so4_a2_anthro-res_surface_mol.nc','w',
                     'NETCDF4_CLASSIC')
    emiss1 = xr.open_dataset('E:/2013_num_so4_a2_anthro-res_surface_mol.nc')