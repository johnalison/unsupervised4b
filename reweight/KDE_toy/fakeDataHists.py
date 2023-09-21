import numpy as np 
import matplotlib.pyplot as plt
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getKDEwJMC, gethist2d, getWeights, getDtoM_loop #getkNNwJCM #, getDtoM_chunk, getDtoMwJCM           
import warnings                                                                                                                                                      
warnings.filterwarnings("ignore", category=RuntimeWarning)
from datetime import datetime


'''- data4bwDtoM	 nSelJets 	 378167.50525067304
- data4bwDtoM	 m4j 	     378167.5052506729
- data4bwDtoM	 leadStM 	 375771.7331247772
- data4bwDtoM	 sublStM 	 378076.50807232986
- data3bwDtoM	 nSelJets 	 8287118.624271587
- data3bwDtoM	 m4j 	     8287121.498107343
- data3bwDtoM	 leadStM 	 8213603.115259366
- data3bwDtoM	 sublStM 	 8283909.236992486'''

def getm4jBinEdges(m,w, bins = 10, m4jMax = 1800):
    m4js = np.where(m<=m4jMax)[0]
    m = m[m4js]; w = w[m4js]
    s = np.argsort(m)
    ws=w[s]; ms=m[s]
    wcum = np.cumsum(ws, dtype=np.float64)
    arr = [np.sum(ws)*i/bins for i in range(1,(bins+1))]
    r = [int(np.round(ms[np.where(wcum<=arr[i])[0][-1]])) for i in range(len(arr))]
    r[-1]=m4jMax; r.insert(0,0)
    bb = [[r[i],r[i+1]] for i in range(len(arr))]
    print(bb)
    np.savetxt('m4jbinbounds2.txt', bb)
    return bb

def dataSel(data, binEdge, dijetLim = None, binEdgeSec = None):
    sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True, False))
    # This is to combine lower and upper sidebands.
    if binEdgeSec is not None:
        sel |= np.asarray(np.where((data['m4j']>binEdgeSec[0]) & (data['m4j']<=binEdgeSec[1]), True, False))

    sel &= np.asarray(np.where((data['leadStM']>0), True,False))
    sel &= np.asarray(np.where((data['sublStM']>0), True,False))
    if dijetLim is not None:
        sel &= np.asarray(np.where((data['leadStM']<=dijetLim), True,False))
        sel &= np.asarray(np.where((data['sublStM']<=dijetLim), True,False))
    return sel


def getFakeQCD(num_events = 1000000, sig_mu = None, sig_std = None, seed = 12345):
    
    def signal(x, sig_mu, sig_std):  
        return np.exp(-((x-sig_mu)/ sig_std)**2) 

    def poly_m4j(x):  
        return np.exp(-((x-500)/250)**2) + 0.4*np.exp(-((x-600)/400)**2)

    def poly_lead(x): 
        return np.exp(-((x-200)/80)**2) + np.exp(-((x-250)/300)**2) + 0.5*np.exp(-((x-100)/300)**2)

    def poly_subl(x):
        return np.exp(-((x-200)/90)**2)  + 0.7*np.exp(-((x-100)/80)**2) + 0.5*np.exp(-((x-100)/300)**2)
    
    def getData( datType = 'm4j', num_events = num_events):
        if datType == 'lead':
            x = np.random.uniform(0, 1000, num_events)
            y = poly_lead(x)
        elif datType == 'subl':
            x = np.random.uniform(0, 1000, num_events)
            y = poly_subl(x)
        elif datType == 'sigm':
            x = np.random.uniform(0, 1800, num_events)
            y = signal(x, 680, 20)
        elif datType == 'sigl':
            x = np.random.uniform(0, 1000, num_events)
            y = signal(x, 91, 10)
        elif datType == 'sigs':
            x = np.random.uniform(0, 1000, num_events)
            y = signal(x, 125, 15)
        else:
            x = np.random.uniform(0, 1800, num_events)
            y = poly_m4j(x)

        return np.random.choice(x, size=num_events, p=y/np.sum(y))
    
    if sig_mu is not None:
        a1, a2, a3 = getData('sigm'), getData('sigl'), getData('sigs')
        b1, b2, b3 = getData('sigm'), getData('sigs'), getData('sigl')
        return np.append(a1,b1), np.append(a2,b2), np.append(a3,b3)
    return getData('m4j'), getData('lead'), getData('subl')

def main(seed = 5, plot = False, signal = False, nBins = 200):

    print("now =", datetime.now())
    m4jBinEdges = np.array([[0, 361], [361, 425], [425, 479], [479, 533], [533, 591], [591, 658], [658, 741], [741, 854], [854, 1044], [1044, 1800]]);

    data4b, data3b = {}, {};

    varRange = [[0,1000],[0,1000]]
    num4b, num3b = int(380000/3), int(8300000/3)
    # num4b, num3b = int(1266), int(27666)
    # num4b, num3b = 100, 2000
    print('Data stats', num4b, num3b)
    np.random.seed(seed)
    data4b['m4j'], data4b['leadStM'], data4b['sublStM'] = getFakeQCD(num_events = num4b)
    data3b['m4j'], data3b['leadStM'], data3b['sublStM'] = getFakeQCD(num_events = num3b)

    m4jBinEdges = np.array(getm4jBinEdges(data3b['m4j'], np.ones(len(data3b['m4j']))))
    
    outputfile = 'hist_noSig_w3to4_toy_seed_'+str(seed)+'_bins_'+str(nBins)
    
    if signal:
        outputfile = 'hist_sig_w3to4_toy_seed_'+str(seed)+'_bins_'+str(nBins)
        
        sel = np.asarray(np.where((data4b['m4j']>m4jBinEdges[7][0]) & (data4b['m4j']<=m4jBinEdges[7][1]), True, False))
        numSig = int(len(data4b['m4j'][sel])*0.05)
        dataSig = {}
        dataSig['m4j'], dataSig['leadStM'], dataSig['sublStM'] = getFakeQCD(num_events = numSig, sig_mu = 680, sig_std = 20, seed = seed)
        numSig = len(dataSig['m4j'])
        num4b +=numSig
        for varType in ['m4j', 'leadStM', 'sublStM']:
            data4b[varType] = np.concatenate((data4b[varType],dataSig[varType]))
        print('updated 4b count',len(data4b[varType]))

    if plot:
        from getQuantilePlots import m4jPlot
            
        data4b['passHLT'] = np.ones(len(data4b['m4j']))
        data3b['passHLT'] = np.ones(len(data3b['m4j']))
        w3to4 = np.load(outputfile+'.npy')
        m4jPlot(m4jBinEdges, 
            data3b = data3b, 
            bgW = w3to4, #np.ones(len(data3b['m4j'])), 
            dataMixed = data4b, 
            dataW = np.ones(len(data4b['m4j'])), 
            figName="test.png", plotAll=False, title='', m4jMax=1800, dijetLim = None, datLabel = None, show=True)

        exit()

    wDtoM4b, wDtoM3b = np.ones(num4b), np.ones(num3b)
    w3to4 = np.zeros(num3b)

    print('Data generated')

    for binNo in range(len(m4jBinEdges)):
        print('SR', binNo, '...', datetime.now())
        if binNo == 0:
            m4jPass_SB3b = dataSel(data3b, m4jBinEdges[binNo+1])
            m4jPass_SB4b = dataSel(data4b, m4jBinEdges[binNo+1])
        elif binNo == len(m4jBinEdges)-1:
            m4jPass_SB3b = dataSel(data3b, m4jBinEdges[binNo-1])
            m4jPass_SB4b = dataSel(data4b, m4jBinEdges[binNo-1])
        else:
            m4jPass_SB3b = dataSel(data3b, m4jBinEdges[binNo-1], binEdgeSec =  m4jBinEdges[binNo+1])
            m4jPass_SB4b = dataSel(data4b, m4jBinEdges[binNo-1], binEdgeSec =  m4jBinEdges[binNo+1])

        leadStM_SB3b = data3b['leadStM'][m4jPass_SB3b]
        sublStM_SB3b = data3b['sublStM'][m4jPass_SB3b]
        # nSelJets_SB3b = data3b['nSelJets'][m4jPass_SB3b]
        wDtoM_SB3b = wDtoM3b[m4jPass_SB3b]
        
        leadStM_SB4b = data4b['leadStM'][m4jPass_SB4b]
        sublStM_SB4b = data4b['sublStM'][m4jPass_SB4b]
        # nSelJets_SB4b = data4b['nSelJets'][m4jPass_SB4b]
        wDtoM_SB4b = wDtoM4b[m4jPass_SB4b]
        
        m4jPass_SR3b = dataSel(data3b, m4jBinEdges[binNo])
        leadStM_SR3b = data3b['leadStM'][m4jPass_SR3b]
        sublStM_SR3b = data3b['sublStM'][m4jPass_SR3b]
        # nSelJets_SR3b = data3b['nSelJets'][m4jPass_SR3b]
        wDtoM_SR3b = wDtoM3b[m4jPass_SR3b]
        
        # num = getDtoM_loop(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, radius=sigma, ttWeight = wDtoM_SB4b)
        # print('num', binNo, '...', datetime.now())
        # den = getDtoM_loop(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, radius=sigma, ttWeight = wDtoM_SB3b)
        # print('den', binNo, '...', datetime.now(), "\n")
        # outputfile = 'w3to4_toy_seed_'+str(seed)+'_radius_'+str(sigma)

        num, binEdgeX, binEdgeY = gethist2d(leadStM_SB4b, sublStM_SB4b, ttWeight = wDtoM_SB4b, nSJtt = None, nBins = nBins, range = varRange)
        den, binEdgeX, binEdgeY = gethist2d(leadStM_SB3b, sublStM_SB3b, ttWeight = wDtoM_SB3b, nSJtt = None, nBins = nBins, range = varRange)

        ratio = np.divide(num, den, out = np.ones(shape=np.shape(den)), where=den!=0)
        w3to4[m4jPass_SR3b] = getWeights(leadStM_SR3b, sublStM_SR3b, ratio, binEdgeX, binEdgeY, w3to4[m4jPass_SR3b])

        # num = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, ttWeight = wDtoM_SB4b, nSJDat = None, nSJtt = None, sigma = sigma, loop = True)
        # den = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, ttWeight = wDtoM_SB3b, nSJDat = None, nSJtt = None, sigma = sigma, loop = True)
        # w3to4[m4jPass_SR3b] = np.divide(num,den)

        # num = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB4b, sigma = sigma)
        # den = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB3b, sigma = sigma)
        # w3to4[m4jPass_SR3b] = np.divide(num,den)

    print("3to4 reweighting complete", datetime.now())
    with uproot.recreate(outputfile+'.root') as outfile:
        outfile["Events"] = {"w3to4": np.array(w3to4)}
        outfile["Events"].show()

    print(w3to4)
    np.save(outputfile, w3to4)
    print(outputfile, 'outputfile saved\n')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument('-4B', '--4b', default='multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0+_newSBDef_2017_multijet.root')
    # parser.add_argument('-3b', '--3B', default='multijet/picoAOD_3b_wJCM_v0_newSBDef_2017_multijet.root')
    parser.add_argument('-s', '--seed', default=2435)
    parser.add_argument('-signal', '--signal', default=False)
    parser.add_argument('-p', '--plot', default=False)
    parser.add_argument('-n', '--nBins', default=15)
    args = parser.parse_args()
    main(int(args.seed), args.plot,bool(args.signal), int(args.nBins))

# import os
# os.system('python3 fakeData.py -s 7 -p True')#

# m4jPlot(m4jBinEdges, 
#             data3b = data3b, 
#             bgW = w3to4, #np.ones(len(data3b['m4j'])), 
#             dataMixed = data4b, 
#             dataW = np.ones(len(data4b['m4j'])), 
#             figName="test.png", plotAll=False, title='seed '+str(seed), m4jMax=1800, dijetLim = None, datLabel = None, show=True)
