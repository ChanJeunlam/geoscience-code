import netCDF4 as nc
import numpy as np

nc_r = nc.Dataset('F:/MEIC/bc_a4_anthro_surface_mol.nc') 
nc_w = nc.Dataset('F:/MEIC/bc_a4_anthro_surface_mol.nc',
                  'w',format = 'NETCDF4')
nc_w.createDimension('date',12) 
nc_w.createDimension('lev',20) 
nc_w.createDimension('ilev',21)
nc_w.createDimension('P0',1)
nc_w.createDimension('lat',192)   
nc_w.createDimension('lon',288)  

##创建变量。参数依次为：‘变量名称’，‘数据类型’，‘基础维度信息’   
nc_w.createVariable('hyai',np.float32,('ilev'))
nc_w.createVariable('hyam',np.float32,('lev'))  
nc_w.createVariable('hybi',np.float32,('ilev'))  
nc_w.createVariable('hybm',np.float32,('lev'))
nc_w.createVariable('ilev',np.float32,('ilev'))
nc_w.createVariable('P0',np.float32,('P0'))    
nc_w.createVariable('lat',np.float32,('lat'))  
nc_w.createVariable('lon',np.float32,('lon'))
nc_w.createVariable('lev',np.float32,('lev'))  
nc_w.createVariable('date',np.int,('date'))
nc_w.createVariable('PS',np.float32,('date','lat','lon')) 
nc_w.createVariable('contvolcano', np.float32, ('date','lev','lat','lon'))

'''
nc_w.createVariable( 'AGR', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'ENE', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'IND', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'RCO', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'SHP', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'SLV', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'TRA', np.float32, ('date','lat','lon'))
nc_w.createVariable( 'WST', np.float32, ('date','lat','lon'))
'''

nc_w.variables['hyai'][:] = nc_r.variables['hyai'][:]
nc_w.variables['hyam'][:] = nc_r.variables['hyam'][:]
nc_w.variables['hybi'][:] = nc_r.variables['hybi'][:]
nc_w.variables['hybm'][:] = nc_r.variables['hybm'][:]
nc_w.variables['ilev'][:] = nc_r.variables['ilev'][:]
nc_w.variables['P0'][:] = nc_r.variables['P0'][:]
nc_w.variables['lat'][:] = nc_r.variables['lat'][:]
nc_w.variables['lon'][:] = nc_r.variables['lon'][:]
nc_w.variables['lev'][:] = nc_r.variables['lev'][:]
nc_w.variables['date'][:] = nc_r.variables['date'][12:24]
nc_w.variables['PS'][:] = nc_r.variables['PS'][12:24]
nc_w.variables['contvolcano'][:] = nc_r.variables['contvolcano'][12:24]

'''
nc_w.variables['AGR'][:] = nc_r.variables['AGR'][2016:2028]
nc_w.variables['ENE'][:] = nc_r.variables['ENE'][2016:2028]
nc_w.variables['IND'][:] = nc_r.variables['IND'][2016:2028]
nc_w.variables['RCO'][:] = nc_r.variables['RCO'][2016:2028]
nc_w.variables['SHP'][:] = nc_r.variables['SHP'][2016:2028]
nc_w.variables['SLV'][:] = nc_r.variables['SLV'][2016:2028]
nc_w.variables['TRA'][:] = nc_r.variables['TRA'][2016:2028]
nc_w.variables['WST'][:] = nc_r.variables['WST'][2016:2028]
'''

nc_w.close()
        
    


