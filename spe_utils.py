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