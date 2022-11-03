import numpy as np
import time
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree

def get_ranges(data, parts = 6):
    count = int(len(data)/parts)
    chunk = []
    for i in range(parts):
        if i == parts-1:
            chunk.append([i*count, len(data)])
        else:
            chunk.append([i*count, (i+1)*count])
    return chunk

def getDtoM(leadData, sublData, leadtt, subltt, radius = 5, ttWeight = None):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    tree = KDTree(ttbar)
    selIndArr = tree.query_radius(data, r = radius)
    if ttWeight is not None:
        return np.asarray([np.sum(ttWeight[i]) for i in selIndArr])
    return np.asarray([len(i) for i in selIndArr])

def getDtoM_chunk(leadData, sublData, leadtt, subltt, radius = 5, ttWeight = None):
    data = np.array([leadData, sublData]).transpose()
    ttbar = np.array([leadtt, subltt]).transpose()
    del leadData, sublData, leadtt, subltt
    tree = KDTree(ttbar)
    bits = get_ranges(data, parts = 60)
    reweight=[]

    if ttWeight is not None:
        for bit in bits:
            selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = 5)
            for i in selIndArr:
                reweight.append(np.sum(ttWeight[i]))
        return np.asarray(reweight)

    for bit in bits:
        selIndArr = tree.query_radius(data[bit[0]:bit[1]], r = 5)
        for i in selIndArr:
            reweight.append(len(i))
    return np.asarray(reweight)
