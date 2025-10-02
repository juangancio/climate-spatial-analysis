import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tqdm, random
from datetime import date, timedelta
import ruptures as rpt
import netCDF4
from statsmodels.tsa.ar_model import AutoReg
from scipy.optimize import curve_fit


def sine(x,p,w,c):
    return p*np.sin(w*x)+c

# Permutation Entropy parameters
L = 4
lag = 1

# Select Dataset and Region

region = 'elnino' # 'elnino' / 'gulf'
d_set = 'ERA5' # 'ERA5' / 'NOAA'

max_CPD_hor = []; max_CPD_ver = []
for n in range(3):

    H_hor = pd.read_csv('surrogates/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'_'+str(n+1)+'.csv',header=None)
    H_ver = pd.read_csv('surrogates/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'_'+str(n+1)+'.csv',header=None)

    #H_hor = pd.read_csv('final_ts/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'.csv',header=None)
    #H_ver = pd.read_csv('final_ts/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'.csv',header=None)
    
    t = np.array(range(len(H_ver)))
    
    '''
    time = np.array(range(len(H_ver)))
    H_hor =np.array(H_hor)[:,0]
    H_ver =np.array(H_ver)[:,0]

    anom = pd.read_csv('final_ts/'+region+'_'+d_set+'_monthly_anomaly.csv',header=None)


    new = H_hor/H_ver

    t1,c =np.polyfit(time,H_hor,1)
    t2,c =np.polyfit(time,H_ver,1)

    H_hor = H_hor-t1*time
    #np.random.shuffle(H_hor)
    #H_hor = H_hor+t1*time
    H_ver = H_hor
    data = H_hor
    model = AutoReg(H_hor, lags=1)
    model_fit = model.fit()
    coef = model_fit.params

    coef_sin = curve_fit(sine,xdata=np.arange(len(data)),ydata=data)[0]

    H_ver = H_ver-t2*time
    #np.random.shuffle(H_ver)
    #H_ver = H_ver+t2*time

    t = np.array(range(len(H_ver)))
    H_ver = sine(t,coef_sin[0],coef_sin[1],coef_sin[2])
    H_hor = [0.67+0.01*np.random.rand(1)]
    for i in range(1,len(t)):
        H_hor.append(H_hor[i-1]*coef[1]+np.random.rand(1)*0+coef[0])

    H_hor = np.array(H_hor)[:,0]
    H_ver = H_ver +np.random.rand(len(t))*np.mean(abs(H_ver-data))
    noise_level = np.mean(abs(H_hor-data))

    for i in range(1,len(t)):
        H_hor[i] = H_hor[i-1]*coef[1]+np.random.rand(1)*noise_level+coef[0]
    '''
    if d_set == 'NOAA':

        # Dataset parameters
        remove_end=1  #some datasets may have not valid values at the end 
        #(usually indicating that the fields have not been computed yet), or it is of insterest that
        #the series ends at a certain point (for example, such as ERA% and NOAA series are of the same length)
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
        remove_end=0   #some datasets may have not valid values at the end 
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


    param = 12
    
    algo = rpt.Pelt(model="rbf").fit(np.array(H_hor))
    result = algo.predict(pen=param)
    results_hor = np.array(result[:-1])

    algo = rpt.Pelt(model="rbf").fit(np.array(H_ver))
    result = algo.predict(pen=param)
    results_ver = np.array(result[:-1])

    '''file_name='final_ts/cpd_pelt_param='+str(param)+'_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_hor.csv'
    np.savetxt(file_name,results_hor, delimiter=",")
    file_name='final_ts/cpd_pelt_param='+str(param)+'_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_ver.csv'
    np.savetxt(file_name,results_ver, delimiter=",")'''

    results_hor = []
    results_ver = []
    param_range = np.linspace(1,16,35)
    time = np.array(time)
    time = t
    for param in tqdm.tqdm(param_range):

        algo = rpt.Pelt(model="rbf").fit(np.array(H_hor))
        #algo = rpt.Window(width=12,model="rbf").fit(np.array(H_hor))
        result = algo.predict(pen=param)

        results_hor.append(np.array(time[result[:-1]]))

        algo = rpt.Pelt(model="rbf").fit(np.array(H_ver))
        #algo = rpt.Window(width=12,model="rbf").fit(np.array(H_ver))
        result = algo.predict(pen=param)
        results_ver.append(np.array(time[result[:-1]]))

    #plt.figure(figsize=(4,16))
    results_hor_array=np.zeros((sst.shape[0],len(param_range)))
    results_ver_array=np.zeros((sst.shape[0],len(param_range)))
    for i in range(len(results_hor)):
        results_hor_array[results_hor[i],i]=param_range[i]
        results_ver_array[results_ver[i],i]=param_range[i]
        if i == 0 and n ==0:
            plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,label='WS',markersize=8)
            plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,label='NS',markersize=8)
        plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,markersize=8)
        plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,markersize=8)

    
    for i in range(sst.shape[0]):
        m = max(results_hor_array[i,:])
        if m>0: max_CPD_hor.append(m)
        m = max(results_ver_array[i,:])
        if m>0: max_CPD_ver.append(m)



#plt.plot(np.linspace(20,25,13),np.linspace(2023,2024,13),'sg',markersize=1)

#plt.yticks(np.array(range(1981,2025)))
plt.grid(); plt.legend(fontsize=12); plt.xlim((min(param_range),max(param_range)))
plt.xlabel('penalty parameter',fontsize=12); plt.ylabel('years',fontsize=12)
plt.xticks(fontsize=12); plt.yticks(fontsize=12)
#plt.savefig('test_elnino.png', format='png', dpi=1200)
#plt.savefig('cpd_region='+region+'_dataset='+d_set+'_lag='+str(lag)+'.png')
plt.show()

label_pos = -18500
label_str = 'Years'
label_pos = -250
label_str = 'Months'

fig = plt.figure(figsize=(6,8))
gs = fig.add_gridspec(2, hspace=.05,height_ratios=[1, 2])
axs = gs.subplots(sharex=True, )
axs[0].plot(time,H_hor)
#axs[0].grid(); axs[0].set_ylabel(r'$H_{WE}$',fontsize=13)
axs[0].text(label_pos,0.7,'(a)',size=13)
axs[0].grid(); axs[0].set_ylabel(r'$AR(1)$',fontsize=13)
ticks = np.array([0.67,0.68,0.69])
axs[0].set_yticks(ticks)
axs[0].set_yticklabels(ticks,fontsize=13)
for i in range(len(results_hor)):
    if i == 0:
        axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sb',alpha=.5,label='WS',markersize=8)
        #axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sr',alpha=.5,label='NS',markersize=8)
        
        #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,label='WS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,label='NS',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,label='WS (sur.)',markersize=8)
        #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,label='NS (sur.)',markersize=8)

    axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sb',alpha=.5,markersize=8)
    #axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sr',alpha=.5,markersize=8)
    

    #plt.plot(param_range[i]*np.ones(len(results_hor[i])),results_hor[i],'sb',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver[i])),results_ver[i],'sr',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_hor_sur[i])),results_hor_sur[i],'sc',alpha=.5,markersize=8)
    #plt.plot(param_range[i]*np.ones(len(results_ver_sur[i])),results_ver_sur[i],'sm',alpha=.5,markersize=8)

#plt.plot(np.linspace(20,25,13),np.linspace(2023,2024,13),'sg',markersize=1)
#plt.text(-5,20000,'(d)',size=13)
#plt.yticks(np.array(range(1981,2025)))
plt.text(label_pos,0.75,'(c)',size=13)
plt.grid(); #plt.legend(fontsize=12); 
plt.ylim((min(param_range),max(param_range)))
plt.ylabel('Penalty parameter',fontsize=13); plt.xlabel(label_str,fontsize=13)
plt.yticks(fontsize=13); plt.xticks(fontsize=13)
plt.gca().invert_yaxis()
#plt.savefig('test_AR.eps', format='eps', dpi=1200, bbox_inches='tight')
#plt.gca().set_rasterized(True)
#plt.gca().set_yticks([1990,2000,2010,2020])
#plt.savefig('cpd_region='+region+'_dataset='+d_set+'_lag='+str(lag)+'_H_hor_surrogate.png', format='png', bbox_inches='tight')
#plt.savefig('cpd_autoregressive.png', format='png', bbox_inches='tight')

plt.show()


fig = plt.figure(figsize=(6,8))
gs = fig.add_gridspec(2, hspace=.05,height_ratios=[1, 2])
axs = gs.subplots(sharex=True, )
axs[0].plot(time,H_ver)
#axs[0].grid(); axs[0].set_ylabel(r'$H_{NS}$',fontsize=13)
axs[0].grid(); axs[0].set_ylabel(r'Noisy sine',fontsize=13)
axs[0].text(label_pos,0.69,'(b)',size=13)
ticks = np.array([0.67,0.68,0.69])
axs[0].set_yticks(ticks)
axs[0].set_yticklabels(ticks,fontsize=13)
for i in range(len(results_hor)):
    if i == 0:
        #axs[1].plot(results_hor[i],param_range[i]*np.ones(len(results_hor[i])),'sb',alpha=.5,label='WS',markersize=8)
        axs[1].plot(results_ver[i],param_range[i]*np.ones(len(results_ver[i])),'sb',alpha=.5,label='NS',markersize=8)
        
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
#plt.text(-5,20000,'(d)',size=13)
#plt.yticks(np.array(range(1981,2025)))
plt.grid(); #plt.legend(fontsize=12); 
plt.ylim((min(param_range),max(param_range)))
plt.ylabel('Penalty parameter',fontsize=13); plt.xlabel(label_str,fontsize=13)
plt.yticks(fontsize=13); plt.xticks(fontsize=13)
plt.gca().invert_yaxis()
plt.text(label_pos,0.75,'(d)',size=13)
#plt.savefig('test_sine.eps', format='eps', dpi=1200, bbox_inches='tight')
#plt.gca().set_rasterized(True)
#plt.gca().set_yticks([1990,2000,2010,2020])
#plt.savefig('cpd_region='+region+'_dataset='+d_set+'_lag='+str(lag)+'_H_ver_surrogate.png', format='png', bbox_inches='tight')
#plt.savefig('cpd_noisy_sine.png', format='png', bbox_inches='tight')

plt.show()

file_name='final_ts/cpd_maxs_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_hor.csv'
#np.savetxt(file_name,max_CPD_hor, delimiter=",")
file_name='final_ts/cpd_max_'+d_set+'_'+region+'_L='+str(L)+'_lag='+str(lag)+'_ver.csv'
#np.savetxt(file_name,max_CPD_ver, delimiter=",")

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