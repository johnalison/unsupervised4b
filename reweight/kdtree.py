import numpy as np
# import time
# import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from sklearn.neighbors import KDTree


def dataSel(data, binEdge, dijetLim = None, binEdgeSec = None):
    sel = np.asarray(np.where(data['passHLT'], True, False))
    sel &= np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True, False))
    # This is to combine lower and upper sidebands.
    if binEdgeSec is not None:
        sel |= np.asarray(np.where((data['m4j']>binEdgeSec[0]) & (data['m4j']<=binEdgeSec[1]), True, False))

    sel &= np.asarray(np.where((data['leadStM']>0), True,False))
    sel &= np.asarray(np.where((data['sublStM']>0), True,False))
    if dijetLim is not None:
        sel &= np.asarray(np.where((data['leadStM']<=dijetLim), True,False))
        sel &= np.asarray(np.where((data['sublStM']<=dijetLim), True,False))
    return sel


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

## xD, yD -> 'data' points that feed KDE, not Data
## xp, yp -> prediction points where KD is estimated
def Gauss2D(xp, yp, xD, yD, sigma, wD = None):
    norm = 1/(2*np.pi*sigma**2)
    gausArr = np.exp(-0.5* (((xp.reshape(-1,1) - xD.reshape(-1))/sigma)**2  + ((yp.reshape(-1,1) - yD.reshape(-1))/sigma)**2  )   )
    if wD is not None:
        return np.sum(gausArr*wD, axis=1)*norm    
    return np.sum(gausArr, axis = 1)*norm

## Calculate the KDE for the given data with sigma = 5 GeV
def getKDE(data, ttbar, ttWeight = None, sigma = 5):
    return Gauss2D(data[:,0], data[:,1], ttbar[:,0], ttbar[:,1], sigma, ttWeight)

def getKDEChunk(data, ttbar, ttWeight = None, sigma = 5, opt=10000):
    opt = int(np.ceil(10000000/len(ttbar[:,0])))
    parts = get_no_of_parts(data[:,0], opt)

    bits = get_ranges(data, parts = parts)
    kdtW=[]

    for bit in bits:
        kdtW.extend(getKDE(data[bit[0]:bit[1]], ttbar, ttWeight, sigma)[:])
    return np.array(kdtW).reshape(-1)


## nSJ = No. of Selected Jets
## JMC = jet multiplicity corrections
def getKDEwJMC(leadData, sublData, leadtt, subltt, ttWeight = None, nSJDat = None, nSJtt = None, sigma = 5):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt

    if nSJDat is not None and ttWeight is not None:
        wArr = np.ones(len(data))
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getKDEChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight[nSJtt==nSJ], sigma=sigma)
        wArr[nSJDat>=9] = getKDEChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight[nSJtt>=9], sigma=sigma )
        return wArr
    elif nSJDat is not None:
        wArr = np.ones(len(data))
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getKDEChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight, sigma=sigma )
        wArr[nSJDat>=9] = getKDEChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight, sigma=sigma )
        return wArr
    else:
        return getKDEChunk(data, ttbar, ttWeight=ttWeight, sigma=sigma)

# get N nearest neighbours of data events and get distance of farthest event
def getNNBall(leadData, sublData, leadtt, subltt, NN = 25, ttWeight = None):
    data = np.array([leadData, sublData]).transpose() # prep the data
    ttbar = np.array([leadtt, subltt]).transpose()
    tree = KDTree(ttbar)
    dist, selIndArr = tree.query(data, k=NN, return_distance=True)
    maxDist = np.array([d[-1] for d in dist])    

    if ttWeight is not None:
        wSum = np.array([np.sum(ttWeight[i]) for i in selIndArr])
        return wSum/np.square(maxDist)  # Return weight density
    countSum = np.array([len(i) for i in selIndArr])
    return countSum/np.square(maxDist)  # Return count density


def getDtoM(leadData, sublData, leadtt, subltt, radius = 5, ttWeight = None):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    tree = KDTree(ttbar)
    selIndArr = tree.query_radius(data, r = radius)
    if ttWeight is not None:
        return np.array([np.sum(ttWeight[i]) for i in selIndArr])
    return np.array([len(i) for i in selIndArr])

def getDtoM_chunk(leadData, sublData, leadtt, subltt, radius = 5, ttWeight = None):
    if len(leadData) <= 10000:
        return getDtoM(leadData, sublData, leadtt, subltt, radius = radius, ttWeight = ttWeight)
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    tree = KDTree(ttbar)
    bits = get_ranges(data, parts = 1000)
    kdtW=[]

    if ttWeight is not None:
        for bit in bits:
            selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = 5)
            for i in selIndArr:
                kdtW.append(np.sum(ttWeight[i]))
        return np.array(kdtW)

    for bit in bits:
        selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = 5)
        for i in selIndArr:
            kdtW.append(len(i))
    return np.array(kdtW)


def getDtoMwJCM(leadData, sublData, leadtt, subltt, radius = 5, ttWeight = None, nSJDat = None, nSJtt = None, parts=100):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    wArr = np.ones(len(data))

    if ttWeight is not None:
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getBallCountChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight[nSJtt==nSJ], parts=parts, radius=radius)
        wArr[nSJDat>=9] = getBallCountChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight[nSJtt>=9], parts=parts, radius=radius)
        return wArr
    for nSJ in range(4,9):
        wArr[nSJDat==nSJ] = getBallCountChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight, parts=parts, radius=radius)
    wArr[nSJDat>=9] = getBallCountChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight, parts=parts, radius=radius)
    return wArr

def getkNNwJCM(leadData, sublData, leadtt, subltt, NN=25, ttWeight = None, nSJDat = None, nSJtt = None, parts=100):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    wArr = np.ones(len(data))

    if ttWeight is not None:
        for nSJ in range(4,9):
            wArr[nSJDat==nSJ] = getkNNDensityChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight[nSJtt==nSJ], parts=parts, NN=NN)
        wArr[nSJDat>=9] = getkNNDensityChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight[nSJtt>=9], parts=parts, NN=3)
        return wArr
    for nSJ in range(4,9):
        wArr[nSJDat==nSJ] = getkNNDensityChunk(data[nSJDat==nSJ], ttbar[nSJtt==nSJ], ttWeight=ttWeight, parts=parts, NN=NN)
    wArr[nSJDat>=9] = getkNNDensityChunk(data[nSJDat>=9], ttbar[nSJtt>=9], ttWeight=ttWeight, parts=parts, NN=3)
    return wArr


def getkNNDensity(data, ttbar, ttWeight=None, NN=25):
    tree = KDTree(ttbar)
    if len(ttbar) < NN:
        NN = int(NN/5)
    dist, selIndArr = tree.query(data, k=NN, return_distance=True)
    maxDist = np.array([d[-1] for d in dist]) 
    if ttWeight is not None:
        wSum = np.array([np.sum(ttWeight[i]) for i in selIndArr])
        return wSum/np.square(maxDist)  # Return weight density
    countSum = np.array([len(i) for i in selIndArr])
    return countSum/np.square(maxDist)  # Return count density


def getkNNDensityChunk(data, ttbar, ttWeight=None, parts=100, NN=25):
    if len(ttbar) < NN and len(ttbar)>5:
        NN = int(NN/5)
    elif len(ttbar) <=5:
        NN = 1
    if len(data) <= 400000:
        return getkNNDensity(data, ttbar, ttWeight, NN)
    tree = KDTree(ttbar)
    bits = get_ranges(data, parts = parts)
    kNND=[]
    if ttWeight is not None:
        for bit in bits:
            dist, selIndArr = tree.query(data[bit[0]:bit[1]], k=NN, return_distance=True)
            maxDist = np.array([d[-1] for d in dist]) 
            for dI, i in enumerate(selIndArr):
                kNND.append(np.sum(ttWeight[i])/np.square(maxDist[dI]))
        return np.array(kNND)
    for bit in bits:
        dist, selIndArr = tree.query(data[bit[0]:bit[1]], k=NN, return_distance=True)
        maxDist = np.array([d[-1] for d in dist]) 
        for dI, i in enumerate(selIndArr):
            kNND.append(len(i)/np.square(maxDist[dI]))
    return np.array(kNND)


def getBallCount(data, ttbar, ttWeight=None, radius=5):
    tree = KDTree(ttbar)
    selIndArr = tree.query_radius(data, r = radius)
    if ttWeight is not None:
        return np.array([np.sum(ttWeight[i]) for i in selIndArr])
    return np.array([len(i) for i in selIndArr])


def getBallCountChunk(data, ttbar, ttWeight=None, parts=100, radius=5):
    if len(data) <= parts*100:
        parts == 1;
        return getBallCount(data, ttbar, ttWeight, radius)
    tree = KDTree(ttbar)
    bits = get_ranges(data, parts = parts)
    kdtW=[]
    if ttWeight is not None:
        for bit in bits:
            selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = radius)
            for i in selIndArr:
                kdtW.append(np.sum(ttWeight[i]))
        return np.array(kdtW)
    for bit in bits:
        selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = radius)
        for i in selIndArr:
            kdtW.append(len(i))
    return np.array(kdtW)