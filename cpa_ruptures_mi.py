import netCDF4
import numpy as np
import math
import matplotlib.pyplot as plt
import random
from datetime import date, timedelta
import scipy
from tqdm import tqdm
import ruptures as rpt
import pandas as pd
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

L = 4
lag = 1

region = 'gulf' # '3' / '4' / '34' / 'gulf'
code = 'horizontal' # 'vertical' / 'horizontal'

#####################################################################################

if region == '34':
    dataset=netCDF4.Dataset(
        "input_data/NOAA_elnino.nc"
        #"/Users/juan/Downloads/oisst_v2_mean_monthly_0.12--150E_-5-5N_1981-2024_0--150E_-4.88-4.88N_1981-2024_160-209.88E_-4.88-4.88N_1981-2024.nc"
        #"/Volumes/T7/climate_data/oisst_v2_anom_monthly_-170-120E_-5-5N_-170--120E_-5-5N.nc"
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
    sst=sst.data[500:1011,0,:,:]

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
    sst2 = sst2.data[500:1011,0,:,:]
    sst = sst.data[500:1011,0,:,:]
    sst = np.concatenate((sst,sst2),axis=2)

elif region == 'gulf':
    dataset=netCDF4.Dataset(
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
time_era = time[500:1026]#[:-2]#
obj = mi_obj(noaa,era,L,lag)
'''obj.code = code
mi = obj.mutal_info()'''
mi_usual = obj.mutal_info_usual()
error = obj.get_error()
pearson_r = obj.get_pearson()

'''mi_usual = mi_usual[310:]
error = error[310:]
pearson_r = pearson_r[310:]'''

time_fit = np.array(range(len(error)))
t1,c =np.polyfit(time_fit,error,1)
error = error - t1*time_fit
plt.plot(error); plt.show()

t2,c =np.polyfit(time_fit,pearson_r,1)
pearson_r = pearson_r - t2*time_fit
plt.plot(pearson_r); plt.show()

t3,c =np.polyfit(time_fit,mi_usual,1)
#mi_usual = mi_usual - t3*time_fit
plt.plot(mi_usual); plt.show()

#a, b = obj.get_corr_length()
#h1, h2 = obj.get_usual_entropies()

#np.savetxt('mi_ts/error_rel_region_'+region+'.csv',mi, delimiter=",")

'''obj.code = 'horizontal'
mi_hor = obj.mutal_info()
obj.code = 'vertical'
mi_ver = obj.mutal_info()'''


mi_hor = pd.read_csv('mi_ts/elnino_anom_MI_hor_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',header=None)
mi_ver = pd.read_csv('mi_ts/elnino_anom_MI_ver_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',header=None)

d_set = 'NOAA'

#mi_hor = pd.read_csv('final_ts/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'.csv',header=None)
#mi_ver = pd.read_csv('final_ts/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'.csv',header=None)



#mi_hor = pd.read_csv('surrogates/elnino_anom_MI_hor_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'_1.csv',header=None)
#mi_ver = pd.read_csv('surrogates/elnino_anom_MI_ver_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'_1.csv',header=None)


mi_hor = np.array(mi_hor); mi_ver = np.array(mi_ver)
mi_hor=mi_hor[:,0]; mi_ver=mi_ver[:,0]; 


time=time_era
time = np.array(time)
t = np.arange(0, time.shape[0])

'''t = t[:480]
time = time[:480]

mi_usual=mi_usual[:480]
error=error[:480]
pearson_r=pearson_r[:480]'''


t = np.array(t)

plt.plot(mi_hor)
t4,c =np.polyfit(t,mi_hor,1)
#mi_hor = mi_hor - t4*t
plt.plot(mi_hor,'r-')
plt.show()

plt.plot(mi_ver)
t5,c =np.polyfit(t,mi_ver,1)
#mi_ver = mi_ver - t5*t
plt.plot(mi_ver,'r-')
plt.show()

t3,c =np.polyfit(t,mi_usual,1)
#mi_usual = mi_usual - t3*t

param = 10
bk_points = 1

#mi_hor = np.column_stack((np.array(mi_hor).reshape(-1, 1), t))
algo = rpt.Pelt(model="rbf").fit(np.array(mi_hor))
result = algo.predict(pen=param)
#result.insert(0,0)
#result[-1] = result[-1] - 1
results_hor = np.array(result[:-1])
print(results_hor)
'''p_value = []; s1 = []; s2 = []; pv1 = []; pv2 = []
for i in range(len(results_hor)-2):
    slope1, intercept, r, p1, se1 = scipy.stats.linregress(t[results_hor[i]:results_hor[i+1]],mi_hor[results_hor[i]:results_hor[i+1],0])
    slope2, intercept, r, p2, se2 = scipy.stats.linregress(t[results_hor[i+1]:results_hor[i+2]],mi_hor[results_hor[i+1]:results_hor[i+2],0])
    z = (slope1-slope2)/np.sqrt(se1**2+se2**2)
    p_value.append(scipy.stats.norm.sf(abs(z))*2)
    s1.append(slope1); s2.append(slope2); pv1.append(p1); pv2.append(p2); 
print(results_hor)
print(time[results_hor])
print(p_value)
print(s1)
print(s2)
print(pv1)
print(pv2)'''

#mi_ver = np.column_stack((np.array(mi_ver).reshape(-1, 1), t))
algo = rpt.Pelt(model="rbf").fit(np.array(mi_ver))
result = algo.predict(pen=param)
#result.insert(0,0)
#result[-1] = result[-1] - 1
results_ver = np.array(result[:-1])
print(results_ver)
'''p_value = []; s1 = []; s2 = []; pv1 = []; pv2 = []
for i in range(len(results_ver)-2):
    slope1, intercept, r, p1, se1 = scipy.stats.linregress(t[results_ver[i]:results_ver[i+1]],mi_ver[results_ver[i]:results_ver[i+1],0])
    slope2, intercept, r, p2, se2 = scipy.stats.linregress(t[results_ver[i+1]:results_ver[i+2]],mi_ver[results_ver[i+1]:results_ver[i+2],0])
    z = (slope1-slope2)/np.sqrt(se1**2+se2**2)
    p_value.append(scipy.stats.norm.sf(abs(z))*2)
    s1.append(slope1); s2.append(slope2); pv1.append(p1); pv2.append(p2); 
print(results_ver)
print(time[results_ver])
print(p_value)
print(s1)
print(s2)
print(pv1)
print(pv2)

results_hor = results_hor[1:-1]
results_ver = results_ver[1:-1]'''

#mi_usual = np.column_stack((np.array(mi_usual).reshape(-1, 1), t))



'''mi_usual = mi_usual[380:]
error = error[380:]
pearson_r = pearson_r[380:]'''
print('Results:')
algo = rpt.Pelt(model="rbf",min_size=24).fit(np.array(mi_usual))
result = algo.predict(pen=param)
results_mi = np.array(result[:-1])
print(results_mi)
#print(time[results_mi])

#error = np.column_stack((np.array(error).reshape(-1, 1), t))
algo = rpt.Pelt(model="rbf",min_size=24).fit(np.array(error))
result = algo.predict(pen=param)
results_error = np.array(result[:-1])
print(results_error)

#pearson_r = np.column_stack((np.array(pearson_r).reshape(-1, 1), t))
algo = rpt.Pelt(model="rbf",min_size=24).fit(np.array(pearson_r))
result = algo.predict(pen=param)
results_pearson = np.array(result[:-1])
print(results_pearson)

file_name='final_ts/cpd_pelt_param='+str(param)+'_SPE_'+region+'_L='+str(L)+'_lag='+str(lag)+'_hor.csv'
#np.savetxt(file_name,results_hor, delimiter=",")
file_name='final_ts/cpd_pelt_param='+str(param)+'_SPE_'+region+'_L='+str(L)+'_lag='+str(lag)+'_ver.csv'
#np.savetxt(file_name,results_ver, delimiter=",")
file_name='final_ts/cpd_pelt_param='+str(param)+'_usual_'+region+'.csv'
#np.savetxt(file_name,results_mi, delimiter=",")
file_name='final_ts/cpd_pelt_param='+str(param)+'_error_'+region+'.csv'
#np.savetxt(file_name,results_error, delimiter=",")
file_name='final_ts/cpd_pelt_param='+str(param)+'_pearson_'+region+'.csv'
#np.savetxt(file_name,results_pearson, delimiter=",")
#print('Results:')
#print(results_mi)

results_hor = []
results_ver = []
results_mi = []
results_error = []
results_pearson = []

param_range = np.linspace(0,200,400)

results_hor_test = np.zeros((len(param_range),len(time)))
results_ver_test = np.zeros((len(param_range),len(time)))
results_mi_test = np.zeros((len(param_range),len(time)))
results_error_test = np.zeros((len(param_range),len(time)))
results_pearson_test = np.zeros((len(param_range),len(time)))


#mi_hor = pd.read_csv('mi_ts/elnino_anom_MI_hor_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',header=None)
#mi_ver = pd.read_csv('mi_ts/elnino_anom_MI_ver_L'+str(L)+'_lag_'+str(lag)+'_region_'+region+'.csv',header=None)

i = 0
for param in tqdm(param_range):
    
    algo = rpt.Pelt(model="rbf").fit(np.array(mi_hor))
    result = algo.predict(pen=param)
    results_hor.append(np.array(time[result[:-1]]))
    results_hor_test[i,result[:-1]] = param

    algo = rpt.Pelt(model="rbf").fit(np.array(mi_ver))
    result = algo.predict(pen=param)
    results_ver.append(np.array(time[result[:-1]]))
    results_ver_test[i,result[:-1]] = param

    algo = rpt.Pelt(model="rbf").fit(np.array(mi_usual))
    result = algo.predict(pen=param)
    results_mi.append(np.array(time[result[:-1]]))
    results_mi_test[i,result[:-1]] = param

    algo = rpt.Pelt(model="rbf").fit(np.array(error))
    result = algo.predict(pen=param)
    results_error.append(np.array(time[result[:-1]]))
    results_error_test[i,result[:-1]] = param

    algo = rpt.Pelt(model="rbf").fit(np.array(pearson_r))
    result = algo.predict(pen=param)
    results_pearson.append(np.array(time[result[:-1]]))
    results_pearson_test[i,result[:-1]] = param

    i =+ 1

m_hor = []; m_ver = []; m_mi = []; m_error = []; m_pearson = []; 
for i in range(len(time)):
    mm_hor = max(results_hor_test[:,i])
    mm_ver = max(results_ver_test[:,i])
    mm_mi = max(results_mi_test[:,i])
    mm_error = max(results_error_test[:,i])
    mm_pearson = max(results_pearson_test[:,i])
    
    if mm_hor>0: m_hor.append(mm_hor)
    if mm_ver>0: m_ver.append(mm_ver)
    if mm_mi>0: m_mi.append(mm_mi)
    if mm_error>0: m_error.append(mm_error)
    if mm_pearson>0: m_pearson.append(mm_pearson)
    

m_hor = np.median(m_hor)
m_ver = np.median(m_ver)
m_mi = np.median(m_mi)
m_error = np.median(m_error)
m_pearson = np.median(m_pearson)


for i in range(len(results_hor)):
    if i == 0:
        plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,label='WE')
        plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,label='NS')
    plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5)
    plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5)

plt.grid(); plt.legend(); plt.xlim((min(param_range),max(param_range)))
plt.xlabel('penalty parameter'); plt.ylabel('years')
#plt.savefig('cpd_region='+region+'_SMI_lag='+str(lag)+'.png')
plt.show()


for i in range(len(results_hor)):
    if i == 0:
        plt.plot(param_range[i]*np.ones(len(results_mi[i])),results_mi[i],'sb',alpha=.5,label='Usual MI')
        plt.plot(param_range[i]*np.ones(len(results_error[i])),results_error[i],'sr',alpha=.5,label='AAD')
        plt.plot(param_range[i]*np.ones(len(results_pearson[i])),results_pearson[i],'sg',alpha=.5,label='pearson r')
        
    plt.plot(param_range[i]*np.ones(len(results_mi[i])),results_mi[i],'sb',alpha=.5)
    plt.plot(param_range[i]*np.ones(len(results_error[i])),results_error[i],'sr',alpha=.5)
    plt.plot(param_range[i]*np.ones(len(results_pearson[i])),results_pearson[i],'sg',alpha=.5)
        
plt.grid(); plt.legend(); plt.xlim((min(param_range),max(param_range)))
plt.xlabel('penalty parameter'); plt.ylabel('years')
#plt.savefig('cpd_region='+region+'_linear_measures.png')
plt.show()



'''plt.plot(time_era,mi)
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
plt.show()'''

'''param = 10

time=time_era

algo = rpt.Pelt(model="rbf").fit(np.array(mi))
result = algo.predict(pen=param)
fig, (ax1, ax2,ax3) = plt.subplots(3)
ax1.plot(time,mi,label='WE symbols')
for i in result[:-1]:
    ax1.plot(time[i],mi[i],'rs',)
plt.ylabel('SPE')
ax1.legend(); ax1.minorticks_on()
ax1.grid(which='both')


ax2.plot(time,error,label='NS symbols')
algo = rpt.Pelt(model="rbf").fit(np.array(error))
result = algo.predict(pen=param)
ax2.legend()
for i in result[:-1]:
    ax2.plot(time[i],error[i],'rs',)
plt.ylabel('SPE')
plt.xlabel('years'); ax2.minorticks_on()
ax2.grid(which='both')

ax3.plot(time,pearson_r,label='NS symbols')
algo = rpt.Pelt(model="rbf").fit(np.array(pearson_r))
result = algo.predict(pen=param)
ax3.legend()
for i in result[:-1]:
    ax3.plot(time[i],pearson_r[i],'rs',)
plt.ylabel('SPE')
plt.xlabel('years'); ax2.minorticks_on()
ax3.grid(which='both')
#plt.savefig('cap_region='+region+'_dataset='+d_set+'_pelt='+str(param)+'_L='+str(L)+'_lag='+str(lag)+'.png', transparent=True); 
plt.show()'''