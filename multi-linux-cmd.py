import os
filename_out = []
cmd =[]
filename_in = os.listdir("/public/home/gaojy/CESM2_MEIC_1850_2017")
for i in range(len(filename_in)):
    filename_out.append('2017_'+filename_in[i])   
    cmd.append('ncks -d time,2016,2027 '+filename_in[i]+' '+filename_out[i])
for i in range(len(filename_in)):
    os.system(cmd[i])

