import numpy as np
# import time
# import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from sklearn.neighbors import KDTree

def get_no_of_parts(data, opt = 10000):
    return(int(np.ceil(len(data)/opt)))

def get_ranges(data, parts = 8):
    count = int(len(data)/parts)
    chunk = []
    for i in range(parts):
        if i == parts-1:
            chunk.append([i*count, len(data)])
        else:
            chunk.append([i*count, (i+1)*count])
    return chunk


######## Needs large memory
## Calculate the KDE for the given data with sigma = 5 GeV
def getKDE(data, ttbar, sigma = 5, ttWeight = None):
    return Gauss2D(data[:,0], data[:,1], ttbar[:,0], ttbar[:,1], sigma, ttWeight)

## xD, yD -> 'data' points that feed KDE, not Data
## xp, yp -> prediction points where KD is estimated
#### Does Broadcasting
def Gauss2D(xp, yp, xD, yD, sigma, wD = None):
    norm = 1/(2*np.pi*sigma**2)
    gausArr = np.exp(-0.5* (((xp.reshape(-1,1) - xD.reshape(-1))/sigma)**2  + ((yp.reshape(-1,1) - yD.reshape(-1))/sigma)**2  )   )
    print(np.shape(gausArr), len(wD))
    if wD is not None:
        return np.sum(gausArr*wD, axis=1)*norm    
    return np.sum(gausArr, axis = 1)*norm

##############################################################################

def Gauss2DLoop(datP, datD, sigma, wD = None):
    if wD is not None:
        return np.sum( np.exp(-0.5*   np.sum(  (( datP - datD) /sigma)**2 , axis = 1) )  *wD)/(2*np.pi*sigma**2)
    return np.sum(   np.exp(-0.5*   np.sum(  (( datP - datD) /sigma)**2 , axis = 1) ) )/(2*np.pi*sigma**2)

def getKDEChunk(data, ttbar, ttWeight = None, sigma = 5, opt=100, loop=True):
    if loop:
        kdtW=[]
        for i in range(int(len(data))):
            kdtW.append(Gauss2DLoop(data[i], ttbar, sigma, ttWeight))
        return np.array(kdtW)
    else:
        return getKDE(data, ttbar, sigma, ttWeight) ######### need large memory
    
    
    # else:
    #     ######## Ignore this ###########
    #     opt = 1
    #     parts= 100
    #     bits = get_ranges(data, parts = parts)
    #     kdtW=[]

    #     for bit in bits:
    #         kdtW.extend(getKDE(data[bit[0]:bit[1]], ttbar, ttWeight, sigma)[:])
    #     return np.array(kdtW).reshape(-1)

## nSJ = No. of Selected Jets
## JMC = jet multiplicity corrections
def getKDEwJMC(leadData, sublData, leadtt, subltt, ttWeight = None, nSJDat = None, nSJtt = None, sigma = 5, loop=True):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt

    if nSJDat is not None and ttWeight is not None:
        wArr = np.ones(len(data))
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getKDEChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight[nSJtt==nSJ], sigma=sigma, loop = loop)
        wArr[nSJDat>=9] = getKDEChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight[nSJtt>=9], sigma=sigma, loop = loop )
        return wArr
    elif nSJDat is not None:
        wArr = np.ones(len(data))
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getKDEChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight, sigma=sigma , loop = loop)
        wArr[nSJDat>=9] = getKDEChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight, sigma=sigma , loop = loop)
        return wArr
    else:
        return getKDEChunk(data, ttbar, ttWeight=ttWeight, sigma=sigma, loop = loop)