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
lag=1


code = 'vertical' # 'vertical' / 'horizontal'

#####################################################################################


dataset=netCDF4.Dataset(
        '/Users/juan/Downloads/data_0.nc'
)

lat = dataset["latitude"][:] #For some reason the data comes N-S inverted
lon = dataset["longitude"][:]
secs = dataset["valid_time"][:]
sst = dataset["t"][:].data
pressure_level = dataset["pressure_level"][:]


base = date(1970, 1, 1)
time=[]
for hour in secs:
    time.append(base + timedelta(seconds=int(hour)))


#sst=space_filter(sst)
l0 = get_anom(sst[:,0,:,:])

mi = []; mi_usual=[]; error=[]; pearson_r=[]; 
for i in range(1,11):
    l1 = get_anom(sst[:,i,:,:])
    #era=sst[:]
    #sst = dataset["sst"][:]
    obj = mi_obj(l0,l1,L,lag)
    obj.code = code
    mi.append(obj.mutal_info())
    mi_usual.append(obj.mutal_info_usual())
    error.append(obj.get_error())
    pearson_r.append(obj.get_pearson())

mi = np.stack([i for i in mi])
mi_usual = np.stack([i for i in mi_usual])
error = np.stack([i for i in error])
pearson_r = np.stack([i for i in pearson_r])

mi = np.swapaxes(mi,0,1)
mi_usual = np.swapaxes(mi_usual,0,1)
error = np.swapaxes(error,0,1)
pearson_r = np.swapaxes(pearson_r,0,1)

#np.savetxt('mi_ts/error_rel_region_'+region+'.csv',mi, delimiter=",")


plt.plot(time,mi)
plt.xlabel('years',size=12); plt.ylabel('SPE-based SMI (a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time,mi_usual)
plt.xlabel('years',size=12); plt.ylabel('usual SMI (a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time,error)
plt.xlabel('years',size=12); plt.ylabel('AAD',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

plt.plot(time,pearson_r)
plt.xlabel('years',size=12); plt.ylabel('r',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
#plt.savefig('H_nino_ERA.png', transparent=True)
plt.show()

'''plt.plot(time,a,label='NOAA'); plt.plot(time_era,b,label='ERA5')
plt.yscale('log')
plt.xlabel('years',size=12); plt.ylabel('correation length(a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
#plt.savefig('corr_length_ver.png', transparent=True)
plt.show()

plt.plot(time,h1,label='NOAA'); plt.plot(time_era,h2,label='ERA5')
plt.yscale('log')
plt.xlabel('years',size=12); plt.ylabel('usual entropy(a.u.)',size=12)
plt.grid(axis='both'); plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
#plt.savefig('corr_length_ver.png', transparent=True)
plt.show()'''

'''if code == 'vertical':  
    np.savetxt('mi_ts/gulf_anom_MI_ver_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',mi, delimiter=",")
elif code == 'horizontal':
    np.savetxt('mi_ts/gulf_anom_MI_hor_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',mi, delimiter=",")
'''