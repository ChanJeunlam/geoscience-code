import xarray as xr

emiss1 = xr.open_dataset('E:/2013_bc_a4_anthro_surface_mol.nc')
emiss2 = xr.open_dataset('E:/2017_bc_a4_anthro_surface_mol.nc')
emiss_mod = xr.open_dataset('E:/ENE_AClean_bc_a4_anthro_surface_mol.nc')
flag_nc = xr.open_dataset('E:/number_29_tag_regions.nc')
for i in range(192):
    for j in range(288):
        if (flag_nc['region_code'].data[i,j]>=61 and 
            flag_nc['region_code'].data[i,j]<=68):
             print(emiss2['ENE'].data[:,i,j])
             print(emiss_mod['ENE'].data[:,i,j])
             print()
                
        if (flag_nc['region_code'].data[i,j]<61 or 
            flag_nc['region_code'].data[i,j]>68):
             print(emiss1['ENE'].data[:,i,j])
             print(emiss_mod['ENE'].data[:,i,j])
             print()