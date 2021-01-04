import pandas
import os
excel_data_df = pandas.read_excel(r'E:\Archive\namelist.xlsx', sheet_name='Sheet1')
namelist_out = excel_data_df['nameID'].tolist()
namelist_in = os.listdir('F:\群相册\群相册')
for i in range(40):
    os.rename('F:/群相册/群相册/'+namelist_in[i], 'F:/群相册/群相册/'+namelist_out[i]+'.png')
