from math import sqrt
from numpy.random import seed
from numpy.random import randn
from numpy import mean
from scipy.stats import t
import netCDF4 as nc
import numpy as np
import plotly.graph_objects as go


# function for calculating the t-test for two dependent samples
def dependent_ttest(data1, data2, alpha):
	# calculate means
	mean1, mean2 = mean(data1), mean(data2)
	# number of paired samples
	n = len(data1)
	# sum squared difference between observations
	d1 = sum([(data1[i]-data2[i])**2 for i in range(n)])
	# sum difference between observations
	d2 = sum([data1[i]-data2[i] for i in range(n)])
	# standard deviation of the difference between means
	sd = sqrt((d1 - (d2**2 / n)) / (n - 1))
	# standard error of the difference between the means
	sed = sd / sqrt(n)
	# calculate the t statistic
	t_stat = (mean1 - mean2) / sed
	# degrees of freedom
	df = n - 1
	# calculate the critical value
	cv = t.ppf(1.0 - alpha, df)
	# calculate the p-value
	p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
	# return everything
	return t_stat, df, cv, p 

modelID = ['CanESM5','CESM2-WACCM','GFDL-ESM4','INM-CM4-8','INM-CM5-0','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']
sspID = ['ssp126','ssp245','ssp370','ssp585']
a = [[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]
a=np.array(a)
b = [[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]
b=np.array(b)
for j in range (0,len(sspID)):
    for i in range (0,len(modelID)):
        nc_obj = nc.Dataset('regrid_emidust_AERmon_'+modelID[i]+'_'+sspID[j]+'_r1i1p1f1_201501-210012.nc')
        emidust = nc_obj.variables['emidust'][:]
        data1 = np.sum(np.sum(np.sum(emidust[0:120, :, :],1),1).reshape(10,12),1)
        data2 = np.sum(np.sum(np.sum(emidust[912:1032, :, :],1),1).reshape(10,12),1)
        
        # calculate the t test
        alpha = 0.05
        t_stat, df, cv, p = dependent_ttest(data1, data2, alpha)
        #print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
        
        # interpret via critical value
        if abs(t_stat) <= cv:
            a[j, i] = 0
        else:
            a[j, i] = 1
        # interpret via p-value
        if p > alpha:
            b[j, i] = 0
        else:
            b[j, i] = 1

fig = go.Figure(data=[go.Table(header=dict(values=[' ','ssp126','ssp245','ssp370','ssp585']),
                 cells=dict(values=[['CanESM5','CESM2-WACCM','GFDL-ESM4','INM-CM4-8','INM-CM5-0','MRI-ESM2-0','NorESM2-LM','NorESM2-MM'],b[0],b[1],b[2],b[3]]))])
fig.show()
