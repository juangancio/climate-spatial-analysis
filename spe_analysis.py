import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
from datetime import date, timedelta
from spe_utils import (
    get_anom,
    h_entropy,
    v_entropy,
    get_mean_anom,
    h_code,
    v_code
)

# Permutation Entropy parameters
L = 4
lag = 8

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
            #"input_data/final_oisst_v2_mean_monthly_-170--120E_5--5N-2.nc"
            "input_data/NOAA_elnino.nc"
                )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            #"input_data/final_oisst_v2_mean_monthly_-67.5--45E_42.5-32.5N.nc"
            "input_data/NOAA_gulf.nc"
                )
    
    elif region == 'sla':    
        dataset = netCDF4.Dataset(
            "input_data/copernicus_sla.nc"
        )

    else:
        raise Exception("Specified region is not valid, try 'elnino' or 'gulf'.")


    time = dataset["time"][:]
    time=time[:-remove_end]
    if region == 'sla':
        sst = dataset["sla"][:]
        time = 1993+time/12
    else:
        time = 1981+9/12+time/12
        sst = dataset["sst"][:]
    sst = sst.data[:-remove_end,::-1,:] #NOAA dataset seems to be inverted NS, this takes dare of that
    #although this does not change the SPE values
    lat = dataset["lat"][:][::-1] #For some reason the data comes N-S inverted
    lon = dataset["lon"][:]

elif d_set == 'ERA5':
    # Dataset parameters
    remove_end=0   #some datasets may have not valid values at the end 
    #(usually indicating that the fields have not been computed yet), or it is of insterest that
    #the series ends at a certain point (for example, such as ERA% and NOAA series are of the same length)
    #this variable takes care of that.

    if region == 'elnino':
        dataset=netCDF4.Dataset(
            #"input_data/final_ERA5_mean_monthly_elnino.nc"
            "input_data/ERA5_elnino.nc"
            )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            #"input_data/final_ERA5_mean_monthly_golfo.nc"
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

#anom = sst[:]
#t_anom= pd.read_csv('final_ts/elnino_'+d_set+'_monthly_anomaly.csv',header=None)


mean_anomaly=[]
std_anomaly=[]
H_hor=[]
H_ver=[]

for i in tqdm.tqdm(range(len(anom))):
    
    #mean_anomaly.append(np.mean(anom[i,:,:]))
    #std_anomaly.append(np.std(anom[i,:,:]))
    H_hor.append(h_entropy(anom[i,:,:],L,lag))
    H_ver.append(v_entropy(anom[i,:,:],L,lag))


t = np.array(range(len(H_ver)))
t1,c =np.polyfit(t,H_hor,1)
t2,c =np.polyfit(t,H_ver,1)
H_hor = H_hor-t1*t
H_ver = H_ver-t2*t


#anom=get_anom(sst[:])
mean_anomaly,std_anomaly = get_mean_anom(anom,lat,lon)

### Save data

#np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anomaly_corrected.csv',mean_anomaly, delimiter=",")
#np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_std_corrected.csv',std_anomaly, delimiter=",")

np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'_detrended_short.csv',H_hor[500:], delimiter=",")
np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'_detrended_short.csv',H_ver[500:], delimiter=",")


### Plotting
plt.plot(time,mean_anomaly)
plt.fill_between(time, 
    np.array(mean_anomaly)-np.array(std_anomaly), 
    np.array(mean_anomaly)+np.array(std_anomaly),
    alpha=0.5)
plt.ylabel('temperature anomaly')
plt.xlabel('years')
plt.show()

plt.plot(time,H_hor,label='WE symbols')
plt.plot(time,H_ver,'r-',label='NS symbols')
plt.ylabel('SPE')
plt.xlabel('years')
plt.legend()
plt.show()
###############################