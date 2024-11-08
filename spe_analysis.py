import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from datetime import date, timedelta
from spe_utils import (
    get_anom,
    h_entropy,
    v_entropy,
    get_mean_anom
)

# Permutation Entropy parameters
L=4
lag=8

# Select Dataset and Region

region = 'elnino' # 'elnino' / 'gulf'
d_set = 'ERA5' # 'ERA5' / 'NOAA'

if d_set == 'NOAA':

    # Dataset parameters
    remove_end=3   #some datasets may have not valid values at the end 
    #(usually indicating that the fields have not been computed yet), or it is of insterest that
    #the series ends at a certain point (for example, such as ERA% and NOAA series are of the same length)
    #this variable takes care of that.
    
    if region == 'elnino':
        dataset=netCDF4.Dataset(
            "input_data/final_oisst_v2_mean_monthly_-170--120E_5--5N-2.nc"
                )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            "input_data/final_oisst_v2_mean_monthly_-67.5--45E_42.5-32.5N.nc"
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

elif d_set == 'ERA5':
    # Dataset parameters
    remove_end=2   #some datasets may have not valid values at the end 
    #(usually indicating that the fields have not been computed yet), or it is of insterest that
    #the series ends at a certain point (for example, such as ERA% and NOAA series are of the same length)
    #this variable takes care of that.

    if region == 'elnino':
        dataset=netCDF4.Dataset(
            "input_data/final_ERA5_mean_monthly_elnino.nc"
            )
    elif region == 'gulf':
        dataset=netCDF4.Dataset(
            "input_data/final_ERA5_mean_monthly_golfo.nc"
                )
    else:
        raise Exception("Specified region is not valid, try 'elnino' or 'gulf'.")


    hours = dataset["time"][:].data

    base = date(1900, 1, 1)
    time=[]
    for hour in hours:
        time.append(base + timedelta(hours=int(hour)))
    time=time[:-remove_end]

    sst = dataset["sst"][:]
    sst=sst.data[:-remove_end,0,:,:]
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

mean_anomaly,std_anomaly = get_mean_anom(anom,lat,lon)

# Save data

np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anomaly_corrected.csv',mean_anomaly, delimiter=",")
np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_std_corrected.csv',std_anomaly, delimiter=",")
np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anom_hor_L'+str(L)+'_lag_'+str(lag)+'.csv',H_hor, delimiter=",")
np.savetxt('final_ts/'+region+'_'+d_set+'_monthly_anom_ver_L'+str(L)+'_lag_'+str(lag)+'.csv',H_ver, delimiter=",")


# Plotting

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