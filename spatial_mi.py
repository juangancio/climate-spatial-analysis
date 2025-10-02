import netCDF4
import numpy as np
import math
import matplotlib.pyplot as plt
import random
from datetime import date, timedelta
import scipy
from tqdm import tqdm
from spe_utils import (
    entropy,
    get_anom,
    h_code,
    v_code
)

random.seed(1)

# CAREFUL!!!
#
# This scripts uses a modified version of the function probabilities, while it usually it
# takes the word length L as a second argument, in this case the number of "bins" ( = L! 
# for ordinal patterns) has to be provided. This was done in order to use the same script 
# to calculate "usual" entropy.
# 
# If the "usual" entropy approach is not used any longer, the script should be corrected
# to be able to work with the same "probabilities" function as the rest of the repository.
# Otherwise, the function here needed should be renamed.

class mi_obj:

    def __init__(self,data1,data2,L,lag):
        
        self.L=L
        self.lag=lag
        #self.lat = int(len(data[0]))
        #self.lon = int(np.size(data[0])/len(data[0]))
        self.data1=data1
        self.data2=data2
        self.code = 'horizontal'
        #self.data=data.reshape(len(data),self.lat*self.lon,1)
        #self.data=self.data[:,:,0]

    def aux_fun(self,i):
       
        if self.code == 'horizontal': 
            var1=h_code(self.data1[i,:,:],self.L,self.lag)
            var2=h_code(self.data2[i,:,:],self.L,self.lag)
        elif self.code == 'vertical':
            var1=v_code(self.data1[i,:,:],self.L,self.lag)
            var2=v_code(self.data2[i,:,:],self.L,self.lag)
        else:
            raise Exception("Specified code is not valid, try 'vertical', or 'horizontal'.")

        h1=entropy(probabilities(var1,math.factorial(self.L)))
        h2=entropy(probabilities(var2,math.factorial(self.L)))
        joint=[var1[i]+(var2[i]-1)*math.factorial(self.L) for i in range(len(var1))]
    
        return h1+h2-entropy(probabilities(joint,math.factorial(self.L)**2))
    
    def aux_fun_usual(self,i):
        
        var1=np.histogram(self.data1[i,:,:].flatten(),bins=math.factorial(self.L),density=True)
        var1 = var1[0]
        h1=scipy.stats.entropy(var1)

        var2 = np.histogram(self.data2[i,:,:].flatten(),bins=math.factorial(self.L),density=True)
        var2 = var2[0]
        h2 = scipy.stats.entropy(var2)

        joint=np.histogram2d(self.data1[i,:,:].flatten(),self.data2[i,:,:].flatten(), bins=math.factorial(self.L),density=True)
        joint=np.ndarray.flatten(joint[0])
        
        return h1+h2-scipy.stats.entropy(joint)
    
    def get_usual_entropies(self):
        h1 = []; h2 = []
        for i in tqdm(range(len(self.data1[:,0,0]))):
            var1=np.histogram(self.data1[i,:,:].flatten(),bins=math.factorial(self.L),density=True)
            var1 = var1[0]
            h1.append(scipy.stats.entropy(var1))

            var2 = np.histogram(self.data2[i,:,:].flatten(),bins=math.factorial(self.L),density=True)
            var2 = var2[0]
            h2.append(scipy.stats.entropy(var2))
        
        return h1, h2
    
    def corr_fun(self,i):
        #var1=v_code(self.data1[i,:,:],self.L,self.lag)
        #var2=v_code(self.data2[i,:,:],self.L,self.lag)
        var1=self.data1[i,:,:].flatten()
        var2=self.data2[i,:,:].flatten()
        return scipy.stats.pearsonr(var1,var2)[0]
    
    def err_fun(self, i):
        var1=self.data1[i,:,:]
        var2=self.data2[i,:,:]
        return np.mean(abs(var1-var2),axis=(0,1))
    
    def mutal_info(self):
        h=[]
        for i in tqdm(range(len(self.data1[:,0,0]))):
            h.append(self.aux_fun(i))
        return h
    def mutal_info_usual(self):
        h=[]
        for i in tqdm(range(len(self.data1[:,0,0]))):
            h.append(self.aux_fun_usual(i))
        return h
    def get_error(self):
        h=[]
        for i in tqdm(range(len(self.data1[:,0,0]))):
            h.append(self.err_fun(i))
        return h
    def get_pearson(self):
        h=[]
        for i in tqdm(range(len(self.data1[:,0,0]))):
            h.append(self.corr_fun(i))
        return h
    
    def get_corr_length(self):
        c1 = []; c2 = [] 
        for i in tqdm(range(len(self.data1[:,0,0]))):
            c_len1 = []; c_len2 = []
            if self.code == 'horizontal':
                for j in range(len(self.data1[0,:,0])):
                    c_len1.append(corr_time(self.data1[i,j,:]))
                    c_len2.append(corr_time(self.data2[i,j,:]))
                c1.append(np.mean(c_len1))
                c2.append(np.mean(c_len2))
            elif self.code == 'vertical':
                for j in range(len(self.data1[0,0,:])):
                    c_len1.append(corr_time(self.data1[i,:,j]))
                    c_len2.append(corr_time(self.data2[i,:,j]))
                c1.append(np.mean(c_len1))
                c2.append(np.mean(c_len2))
            else:
                raise Exception("Specified code is not valid, try 'vertical', or 'horizontal'.")
        return c1, c2



def probabilities(code,bins):# CAREFUL!!! This is not the usual function!
    get_indexes = lambda x, xs: [k for (y, k) in zip(xs, range(len(xs))) if x == y]
    probs=[]
    for i in range(1,bins+1):
 
        probs=probs + [len(get_indexes(i,code))/len(code)
                   ]
        
    return probs


def usual_entropy(data,b):
    pk=np.histogram(data,bins=b,density=True)
    return scipy.stats.entropy(pk[0])/np.log(b)

def corr_time(x):
    x = x-np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[result.size//2:]

    result = result/np.mean(x**2)
    result = result**2
    return scipy.integrate.trapezoid(result,dx=1)

#####################################################################################

L=4
lag=8

region = 'gulf' # '3' / '4' / '34' / 'gulf'
code = 'vertical' # 'vertical' / 'horizontal'

#####################################################################################

if region == '34':
    dataset=netCDF4.Dataset(
        #"input_data/final_oisst_v2_mean_monthly_-170--120E_5--5N-2.nc"
        #"/Users/juan/Downloads/oisst_v2_mean_monthly_0.12--150E_-5-5N_1981-2024_0--150E_-4.88-4.88N_1981-2024_160-209.88E_-4.88-4.88N_1981-2024.nc"
        #"/Volumes/T7/climate_data/oisst_v2_anom_monthly_-170-120E_-5-5N_-170--120E_-5-5N.nc"
        "input_data/NOAA_elnino.nc"
        )
elif region == '3':
    dataset=netCDF4.Dataset(
        "input_data/oisst_v2_mean_monthly_-150--90E_-5-5N.nc"
    )
elif region == '4':
    dataset=netCDF4.Dataset(
        "input_data/oisst_v2_mean_monthly_0.12--150E_-5-5N_1981-2024_0--150E_-4.88-4.88N_1981-2024_160-209.88E_-4.88-4.88N_1981-2024.nc"
    )
elif region == 'gulf':
    dataset=netCDF4.Dataset(
        #"input_data/final_oisst_v2_mean_monthly_-67.5--45E_42.5-32.5N.nc"
        "input_data/NOAA_gulf.nc"
    )  
else:
    raise Exception("Specified region is not valid, try '3', '4', '34', or 'gulf.")

lat_sat = dataset["lat"][:][::-1] #For some reason the data comes N-S inverted
lon_sat = dataset["lon"][:]
time_sat = dataset["time"][:]
sst = dataset["sst"][:]
anom=get_anom(sst[:])
#noaa=sst[:]
#anom = dataset["sla"][:]
anom = anom[:,::-1,:] #For some reason the data comes N-S inverted
noaa=anom.data[:,:,:]
time_sat=time_sat[:526]
noaa=noaa[:526,:,:]

print('NOAA data loaded')

if region == '34':
    dataset=netCDF4.Dataset(
        #"input_data/final_ERA5_mean_monthly_elnino.nc"
        "input_data/ERA5_elnino.nc"
        #"/Users/juan/Downloads/adaptor.mars.internal-1724666270.488026-15195-15-b85dba3b-0e86-4e89-b008-721e5fb0829d.nc"
        #"/Volumes/T7/climate_data/adaptor.mars.internal-1716106111.9322994-1587-14-4a59a5e2-b6ae-4d38-a377-800bb134c4e8.nc"
        )

    lon_era = dataset["longitude"][:]
    sst = dataset["sst"][:]
    sst=sst.data[500:1026,:,:]

elif region == '3':
    dataset=netCDF4.Dataset(
        "input_data/ERA5_mean_monthly_elnino3.nc"
    )
    
    lon_era = dataset["longitude"][:]
    sst = dataset["sst"][:]
    sst=sst.data[500:1026,:,:]

elif region == '4':
    dataset=netCDF4.Dataset(
        "input_data/adaptor.mars.internal-1724686270.1487467-12125-1-cca537d8-9827-4324-92ce-3f95fb19035c.nc"
    )
    dataset2=netCDF4.Dataset(
        "input_data/adaptor.mars.internal-1724745632.9268458-3371-18-14c0ba75-67c9-4043-ac8b-58a389fae2e9.nc"
    )
    lon_era=np.concatenate((dataset["longitude"][:]-360,dataset2["longitude"][:]))
    sst = dataset["sst"][:]
    sst2 = dataset2["sst"][:]
    sst2 = sst2.data[500:1026,0,:,:]
    sst = sst.data[500:1026,0,:,:]
    sst = np.concatenate((sst,sst2),axis=2)

elif region == 'gulf':
    dataset=netCDF4.Dataset(
        #"input_data/final_ERA5_mean_monthly_golfo.nc"
        "input_data/ERA5_gulf.nc"
    )
    lon_era = dataset["longitude"][:]
    sst = dataset["sst"][:]
    sst=sst.data[500:1026,:,:]

else:
    raise Exception("Specified region is not valid, try '3', '4', '34', or 'gulf'.")

hours = dataset["valid_time"][:].data
lat_era = dataset["latitude"][:]

base = date(1970, 1, 1)
time=[]
for hour in hours:
    time.append(base + timedelta(seconds=int(hour)))

print('ERA5 data loaded')


#sst=space_filter(sst)
era = get_anom(sst[:])
#era=sst[:]
#sst = dataset["sst"][:]
time_era = time[500:1026]
obj = mi_obj(noaa,era,L,lag)
obj.code = code
mi = obj.mutal_info()
mi_usual = obj.mutal_info_usual()
error = obj.get_error()
pearson_r = obj.get_pearson()
a, b = obj.get_corr_length()
h1, h2 = obj.get_usual_entropies()

#np.savetxt('mi_ts/usual_mi_region_'+region+'.csv',mi_usual, delimiter=",")
#np.savetxt('mi_ts/error_region_'+region+'.csv',error, delimiter=",")
#np.savetxt('mi_ts/pearson_region_'+region+'.csv',pearson_r, delimiter=",")

plt.plot(time_era,mi)
plt.xlabel('years',size=12); plt.ylabel('SPE-based SMI (a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time_era,mi_usual)
plt.xlabel('years',size=12); plt.ylabel('usual SMI (a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time_era,error)
plt.xlabel('years',size=12); plt.ylabel('AAD',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time_era,pearson_r)
plt.xlabel('years',size=12); plt.ylabel('r',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time_era,a,label='NOAA'); plt.plot(time_era,b,label='ERA5')
plt.yscale('log')
plt.xlabel('years',size=12); plt.ylabel('correation length(a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
#plt.savefig('corr_length_ver.png', transparent=True)
plt.show()

plt.plot(time_era,h1,label='NOAA'); plt.plot(time_era,h2,label='ERA5')
plt.yscale('log')
plt.xlabel('years',size=12); plt.ylabel('usual entropy(a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
#plt.savefig('corr_length_ver.png', transparent=True)
plt.show()

if code == 'vertical':  
    np.savetxt('mi_ts/elnino_anom_MI_ver_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',mi, delimiter=",")
elif code == 'horizontal':
    np.savetxt('mi_ts/elnino_anom_MI_hor_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',mi, delimiter=",")
