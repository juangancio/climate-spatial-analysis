import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import tqdm
from datetime import date, timedelta
import ruptures as rpt
from spe_utils import (
    get_anom,
    h_entropy,
    v_entropy,
    get_mean_anom
)

# Permutation Entropy parameters
L = 4
lag = 1

# Select Dataset and Region

region = 'gulf' # 'elnino' / 'gulf'
d_set = 'ERA5' # 'ERA5' / 'NOAA'

if d_set == 'NOAA':

    # Dataset parameters
    remove_end=1   #some datasets may have not valid values at the end 
    #(usually indicating that the fields have not been computed yet), or it is of insterest that
    #the series ends at a certain point (for example, such as ERA5 and NOAA series are of the same length)
    #this variable takes care of that.
    
    if region == 'elnino':
        dataset=netCDF4.Dataset(
            "input_data/NOAA_elnino.nc"
                )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            "input_data/NOAA_gulf.nc"
                )
    else:
        raise Exception("Specified region is not valid, try 'elnino' or 'gulf'.")


    time = dataset["time"][:]
    time=time[:-remove_end]
    sst = dataset["sst"][:]
    sst = sst.data[:-remove_end,::-1,:] #NOAA dataset seems to be inverted NS, this takes dare of that
    #although this does not change the SPE values
    lat = dataset["lat"][:][::-1] #For some reason the data comes N-S inverted
    lon = dataset["lon"][:]

    time = 1981+9/12+time/12

elif d_set == 'ERA5':
    # Dataset parameters
    remove_end=2   #some datasets may have not valid values at the end 
    #(usually indicating that the fields have not been computed yet), or it is of insterest that
    #the series ends at a certain point (for example, such as ERA% and NOAA series are of the same length)
    #this variable takes care of that.

    if region == 'elnino':
        dataset=netCDF4.Dataset(
            "input_data/ERA5_elnino.nc"
            )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            "input_data/ERA5_gulf.nc"
                )
    else:
        raise Exception("Specified region is not valid, try 'elnino' or 'gulf'.")


    hours = dataset["valid_time"][:].data

    base = date(1970, 1, 1)
    time=[]
    for hour in hours:
        time.append(base + timedelta(seconds=int(hour)))
    time=time[:]

    sst = dataset["sst"][:]
    sst=sst.data[:,:,:]
    lon = dataset["longitude"][:]
    lat = dataset["latitude"][:]

else:
    raise Exception("Specified dataset is not valid, try 'ERA5' or 'NOAA'.")


anom=get_anom(sst[:])

mean_anomaly=[]
std_anomaly=[]
H_hor=[]
H_ver=[]

for i in tqdm.tqdm(range(len(anom))):
    
    #mean_anomaly.append(np.mean(anom[i,:,:]))
    #std_anomaly.append(np.std(anom[i,:,:]))
    H_hor.append(h_entropy(anom[i,:,:],L,lag))
    H_ver.append(v_entropy(anom[i,:,:],L,lag))

'''H_hor = pd.read_csv('surrogates/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'_3.csv',header=None)
H_ver = pd.read_csv('surrogates/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'_3.csv',header=None)

H_hor = np.array(H_hor)[:,0]; H_ver = np.array(H_ver)[:,0]; '''

'''H_hor = H_hor[500:]
H_ver = H_ver[500:]
time = time[500:]'''

t = np.array(range(len(H_ver)))
t1,c =np.polyfit(t,H_hor,1)
t2,c =np.polyfit(t,H_ver,1)
H_hor = H_hor-t1*t
H_ver = H_ver-t2*t

H_hor = np.array(H_hor)
H_ver = np.array(H_ver)

mean_anomaly,std_anomaly = get_mean_anom(anom,lat,lon)

file_name='final_ts/cpd_maxs_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_hor.csv'
maxs_hor = pd.read_csv(file_name,header=None)
file_name='final_ts/cpd_max_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_ver.csv'
maxs_ver = pd.read_csv(file_name,header=None)
#print(len(maxs_hor)); print(len(maxs_ver))

thrs_hor1 = np.percentile(maxs_hor,99.5)
thrs_ver1 = np.percentile(maxs_ver,99.5)

thrs_hor2 = np.percentile(maxs_hor,95)
thrs_ver2 = np.percentile(maxs_ver,95)

param = thrs_hor1
algo = rpt.Pelt(model="rbf").fit(np.array(H_hor))
result = algo.predict(pen=param)
results_hor = np.array(result[:-1])

param = thrs_ver1
algo = rpt.Pelt(model="rbf").fit(np.array(H_ver))
result = algo.predict(pen=param)
results_ver = np.array(result[:-1])

file_name='final_ts/cpd_pelt_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_hor.csv'
np.savetxt(file_name,results_hor, delimiter=",")
file_name='final_ts/cpd_pelt_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_ver.csv'
np.savetxt(file_name,results_ver, delimiter=",")
print(results_hor)
print(results_ver)

results_hor = []
results_ver = []
results_hor_sur = []
results_ver_sur = []

param_range = np.linspace(0,50,100)
time = np.array(time)

#H_hor = np.column_stack((np.array(H_hor).reshape(-1, 1), t))
#H_ver = np.column_stack((np.array(H_ver).reshape(-1, 1), t))


results_hor_test = np.zeros((len(param_range),len(time)))
results_ver_test = np.zeros((len(param_range),len(time)))

i = 0
for param in tqdm.tqdm(param_range):
    
    algo = rpt.Pelt(model="rbf").fit(np.array(H_hor))
    #algo = rpt.Dynp(model="clinear").fit(H_hor)
    result = algo.predict(pen=param)
    results_hor.append(np.array(time[result[:-1]]))
    results_hor_test[i,result[:-1]] = param
    
    algo = rpt.Pelt(model="rbf").fit(np.array(H_ver))
    #algo = rpt.Dynp( model="clinear").fit(H_ver)
    result = algo.predict(pen=param)
    results_ver.append(np.array(time[result[:-1]]))
    results_ver_test[i,result[:-1]] = param

    i =+ 1
    #algo = rpt.Pelt(model="rbf").fit(np.array(H_hor_sur))
    #result = algo.predict(pen=param)
    #results_hor_sur.append(np.array(time[result[:-1]]))

    #algo = rpt.Pelt(model="rbf").fit(np.array(H_ver_sur))
    #result = algo.predict(pen=param)
    #results_ver_sur.append(np.array(time[result[:-1]]))

m_hor = []; m_ver = []
for i in range(len(time)):
    mm_hor = max(results_hor_test[:,i])
    mm_ver = max(results_ver_test[:,i])
    if mm_hor>0: m_hor.append(mm_hor)
    if mm_ver>0: m_ver.append(mm_ver)

m_hor = np.median(m_hor)
m_ver = np.median(m_ver)

label_pos = -18500
#label_str = 'Years'
label_pos = 0#-250
#label_str = 'Months'

#plt.figure(figsize=(4,16))
fig = plt.figure(figsize=(6,8))
gs = fig.add_gridspec(2, hspace=.05,height_ratios=[1, 2])
axs = gs.subplots(sharex=True, )
axs[0].plot(time,H_hor)
plt.figtext(label_pos,0.87,'(a)',size=13)
axs[0].grid(); #axs[0].set_ylabel(r'$H_{WE}$',fontsize=14)
axs[0].set_ylabel(r'$H_{WE} (detrended)$',fontsize=13)
ticks = np.array([0.6,0.65,0.7,0.75])
axs[0].set_yticks(ticks)
axs[0].set_yticklabels(ticks,fontsize=13)
for i in range(len(results_hor)):
    if i == 0:
        axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sc',alpha=.5,markersize=8)
        #axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sr',alpha=.5,label='NS',markersize=8)
        
        #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,label='WS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,label='NS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,label='WS (sur.)',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,label='NS (sur.)',markersize=8)

    axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sc',alpha=.5,markersize=8)
    #axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sr',alpha=.5,markersize=8)
    

    #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,markersize=8)

#plt.plot(np.linspace(20,25,13),np.linspace(2023,2024,13),'sg',markersize=1)
#plt.text(-5,20000,'(d)',size=14)
#plt.yticks(np.array(range(1981,2025)))
plt.plot([time[-1],time[0]],[thrs_hor1,thrs_hor1],'k-',label=r'99.5th percentile of $P^*$ from surrogates')
plt.plot([time[-1],time[0]],[m_hor,m_hor],'k--',label=r'$\overline{P^{*}}$ of the signal')
plt.grid(); #plt.legend(fontsize=12); 
plt.figtext(label_pos,0.6,'(c)',size=13)

'''plt.figtext(.75,.2,r'$P^*(2007)$',size=13)
plt.figtext(.5,.47,r'$P^*(1982)$',size=13)
plt.figtext(.65,.35,r'$R=51.7$',size=13, rotation='vertical')
plt.figtext(.45,.5,r'$R=13.7$',size=13, rotation='vertical')
plt.annotate('', xy=(12000, 1), xytext=(12000, 32),
            arrowprops=dict(arrowstyle='<->', color='red'))
plt.annotate('', xy=(3500, 1), xytext=(3500, 11),
            arrowprops=dict(arrowstyle='<->', color='red'))
plt.plot(time[875], max(results_hor_test[:,875]),'sr',markersize=8)
plt.plot(time[810], max(results_hor_test[:,810]),'sr',markersize=8)'''

plt.figtext(.6,.26,r'$P^*(2007)$',size=13)
plt.figtext(.76,.26,r'$P^*(2013)$',size=13)
plt.figtext(.66,.4,r'$R=61.0$',size=13, rotation='vertical')
plt.figtext(.8,.4,r'$R=61.0$',size=13, rotation='vertical')
plt.annotate('', xy=(12500, 1), xytext=(12500, 32),
            arrowprops=dict(arrowstyle='<->', color='red'))
plt.annotate('', xy=(17000, 1), xytext=(17000, 32),
            arrowprops=dict(arrowstyle='<->', color='red'))
plt.plot(time[875], max(results_hor_test[:,875]),'sr',markersize=8)
plt.plot(time[810], max(results_hor_test[:,810]),'sr',markersize=8)

plt.ylim((min(param_range),max(param_range)))
plt.ylabel('Penalty parameter',fontsize=13); plt.xlabel('Years',fontsize=13)
plt.yticks(fontsize=14); plt.xticks(fontsize=13)
plt.gca().invert_yaxis()
plt.legend(loc='lower left',fontsize=11)
plt.savefig('test_iAAFT2_v2.eps', format='eps', bbox_inches='tight')
#plt.gca().set_rasterized(True)
#plt.gca().set_yticks([1990,2000,2010,2020])
#plt.savefig('cpd_region='+region+'_dataset='+d_set+'_lag='+str(lag)+'_H_hor_short.png', format='png', bbox_inches='tight')

plt.show()


fig = plt.figure(figsize=(6,8))
gs = fig.add_gridspec(2, hspace=.05,height_ratios=[1, 2])
axs = gs.subplots(sharex=True, )
axs[0].plot(time,H_ver)
plt.figtext(label_pos,0.87,'(b)',size=13)
axs[0].grid(); #axs[0].set_ylabel(r'$H_{NS}$',fontsize=14)
axs[0].set_ylabel(r'$H_{NS}$ (detrended)',fontsize=14)
#axs[0].set_yticks([0.5,0.55,0.6,0.65])
#axs[0].set_yticklabels([0.5,0.55,0.6,0.65],fontsize=14)
for i in range(len(results_hor)):
    if i == 0:
        #axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sb',alpha=.5,label='WS',markersize=8)
        axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sb',alpha=.5,markersize=8)
        
        #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,label='WS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,label='NS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,label='WS (sur.)',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,label='NS (sur.)',markersize=8)

    #axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sb',alpha=.5,markersize=8)
    axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sb',alpha=.5,markersize=8)
    

    #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,markersize=8)

#plt.plot(np.linspace(20,25,13),np.linspace(2023,2024,13),'sg',markersize=1)
#plt.text(-5,20000,'(d)',size=14)
#plt.yticks(np.array(range(1981,2025)))
plt.plot([time[-1],time[0]],[thrs_ver1,thrs_ver1],'r-', label='99.5th percentile')
#plt.plot([time[-1],time[0]],[thrs_ver2,thrs_ver2],'g-', label='95th percentile')
plt.grid(); #plt.legend(fontsize=12); 
plt.figtext(label_pos,0.6,'(d)',size=13)
plt.ylim((min(param_range),max(param_range)))
plt.ylabel('Penalty parameter',fontsize=14); plt.xlabel('Years',fontsize=14)
plt.yticks(fontsize=14); plt.xticks(fontsize=14)
plt.gca().invert_yaxis()
plt.legend(loc='lower left',fontsize=13)
#plt.savefig('test_elnino.png', format='png', dpi=1200)
#plt.gca().set_rasterized(True)
#plt.gca().set_yticks([1990,2000,2010,2020])
#plt.savefig('cpd_region='+region+'_dataset='+d_set+'_lag='+str(lag)+'_H_ver_short.png', format='png', bbox_inches='tight')

plt.show()

'''fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(time,H_hor,label='WE symbols')
for i in result[:-1]:
    ax1.plot(time[i],H_hor[i],'rs',)
plt.ylabel('SPE')
ax1.legend(); ax1.minorticks_on()
ax1.grid(which='both')


ax2.plot(time,H_ver,label='NS symbols')
algo = rpt.Pelt(model="rbf").fit(np.array(H_ver))
result = algo.predict(pen=param)
ax2.legend()
for i in result[:-1]:
    ax2.plot(time[i],H_ver[i],'rs',)
plt.ylabel('SPE')
plt.xlabel('years'); ax2.minorticks_on()
ax2.grid(which='both')
#plt.savefig('cap_region='+region+'_dataset='+d_set+'_pelt='+str(param)+'_L='+str(L)+'_lag='+str(lag)+'.png', transparent=True); 
plt.show()'''
###############################