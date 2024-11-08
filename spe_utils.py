import numpy as np
import math


def get_anom(sst):
    lat = int(len(sst[0]))
    lon = int(np.size(sst[0])/len(sst[0]))
    ave = np.array([[[0] * lon] * lat]*12,dtype='float64')

    for i in range(lat):
        for j in range(lon):
            for month  in range(12):
                ave[month,i,j]=np.mean([sst[t,i,j] for t in range(month,len(sst),12)])
    
    anom=sst[:,:,:]
    month=0
    for t in range(len(sst)):
        for i in range(lat):
            for j in range (lon):
                anom[t,i,j]=anom[t,i,j]-ave[month,i,j]
        month+=1
        if month==12:
            month=0
    return anom

def get_mean_anom(anom,lat,lon):
    anom = np.array(anom)
    grids = np.meshgrid(lon,lat)
    grid_lat = np.array(grids[0]); grid_lon = np.array(grids[1]) 

    areas = grid_area(grid_lat,grid_lon)
    ave = []; std = []
    for i in range(anom.shape[0]):
        #ave.append(np.sum(anom[i,:,:]*areas))
        mean = np.average(anom[i,:,:],weights=areas)
        var = np.average((anom[i,:,:]-mean)**2, weights=areas)
        ave.append(mean); std.append(np.sqrt(var))
    return (np.array(ave),std)#/np.sum(areas)

def earth_radius(lat,lon): #Eartch radius in meters
    a = 6378137 # equatorial radius
    b = 6356752 # polar radius
    r = (((a**2*np.cos(np.radians(lat)))**2 + (b**2*np.sin(np.radians(lat)))**2) /
        ((a*np.cos(np.radians(lat)))**2 + (b*np.sin(np.radians(lat)))**2))**(1/2)
    return r

def grid_area(lat,lon):

    R = earth_radius(lat,lon) # Earth radius in meters

    ## Determine grid sizes dlat and dlon: 
    [dlat1,dlat2] = np.gradient(lat); 
    [dlon1,dlon2] = np.gradient(lon); 
    # We don't know if lat and lon were created by [lat,lon] = meshgrid(lat_array,lon_array) or [lon,lat] = meshgrid(lon_array,lat_array) 
    # but we can find out: 
    if np.array_equal(dlat1,np.zeros(lat.shape)):
        dlat = dlat2; 
        dlon = dlon1; 
        if not(np.array_equal(dlon2, np.zeros(lon.shape))):
            raise Exception('Error: lat and lon must be monotonic grids, as if created by meshgrid.') 
    else:
        dlat = dlat1; 
        dlon = dlon2; 
        if not(np.array_equal(dlon1,dlat2)) or not(np.array_equal(dlat2,np.zeros(lon.shape))) or not(np.array_equal(dlon1,np.zeros(lon.shape))):
            raise Exception('Error: lat and lon must be monotonic grids, as if created by meshgrid.') 
    
    ## Calculate area based on dlat and dlon: 
    dy = dlat*R*np.pi/180
    dx = (dlon/180)*np.pi*R*np.cos(np.radians(lat)); 
    return abs(dx*dy)


def h_entropy(data,L,lag):
    
    code=[]
    for i in range(len(data[:,0])):

        code.append(perm_indices(data[i,:],L,lag))
    
    return entropy(probabilities(np.array(code).reshape(1,np.size(code))[0],L))/np.log(math.factorial(L))

def v_entropy(data,L,lag):
    
    code=[]
    for i in range(len(data[0,:])):

        code.append(perm_indices(data[:,i],L,lag))
    
    return entropy(probabilities(np.array(code).reshape(1,np.size(code))[0],L))/np.log(math.factorial(L))

def perm_indices(ts, wl, lag):
    m = len(ts) - (wl - 1)*lag
    indcs = np.zeros(m, dtype=int)
    for i in range(1,wl):
        st = ts[(i - 1)*lag : m + ((i - 1)*lag)]
        for j in range(i,wl):
            zipped=zip(st,ts[j*lag : m+j*lag])
            indcs += [x > y for (x, y) in zipped]

        indcs*= wl - i
    return indcs + 1

def h_code(data,L,lag):
    code=[]
    for i in range(len(data[:,0])):
        code.append(perm_indices(data[i,:],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return np.array(code).reshape(1,np.size(code))[0]

def v_code(data,L,lag):
    code=[]
    for i in range(len(data[0,:])):
        code.append(perm_indices(data[:,i],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return np.array(code).reshape(1,np.size(code))[0]




def entropy(probs):
    h=0
    for i in range(len(probs)):
        if probs[i]==0:
            continue
        else:
            h=h-probs[i]*np.log(probs[i])
    return h

def probabilities(code,L):
    get_indexes = lambda x, xs: [k for (y, k) in zip(xs, range(len(xs))) if x == y]
    probs=[]
    for i in range(1,math.factorial(L)+1):
 
        probs=probs + [len(get_indexes(i,code))/len(code)
                   ]
        
    return probs