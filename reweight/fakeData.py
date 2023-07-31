import numpy as np 
import matplotlib.pyplot as plt
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getKDEwJMC #getkNNwJCM #, getDtoM_chunk, getDtoMwJCM               
import awkward as ak
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


def getFakeQCD(num_events = 1000000, seed = 5):
    
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
        else:
            x = np.random.uniform(0, 1800, num_events)
            y = poly_m4j(x)

        return np.random.choice(x, size=num_events, p=y/np.sum(y))
    
    np.random.seed(seed)
    return getData('m4j'), getData('lead'), getData('subl')
    

def main(seed, plot = False):

    print("now =", datetime.now())
    m4jBinEdges = np.array([[0, 361], [361, 425], [425, 479], [479, 533], [533, 591], [591, 658], [658, 741], [741, 854], [854, 1044], [1044, 1800]]);

    seed = int(seed)
    sigma = 5; 
    outputfile = 'w3to4_toy_seed_'+str(seed)+'_sigma_'+str(sigma)
    data4b, data3b = {}, {};

    num4b, num3b = int(380000/3), int(8300000/3)
    # num4b, num3b = 10, 200
    print('Data stats', num4b, num3b)
    data4b['m4j'], data4b['leadStM'], data4b['sublStM'] = getFakeQCD(num_events = num4b, seed = seed)
    data3b['m4j'], data3b['leadStM'], data3b['sublStM'] = getFakeQCD(num_events = num3b, seed = seed)


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

        num = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, ttWeight = wDtoM_SB4b, nSJDat = None, nSJtt = None, sigma = sigma, loop = True)
        print('num', binNo, '...', datetime.now())
        den = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, ttWeight = wDtoM_SB3b, nSJDat = None, nSJtt = None, sigma = sigma, loop = True)
        print('den', binNo, '...', datetime.now(), "\n")
        # num = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB4b, sigma = sigma)
        # den = getKDEwJMC(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB3b, sigma = sigma)
        w3to4[m4jPass_SR3b] = np.divide(num,den)

    print("3to4 reweighting complete", datetime.now())
    np.save(outputfile, w3to4)
    print(outputfile, 'outputfile saved\n')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument('-4B', '--4b', default='multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0+_newSBDef_2017_multijet.root')
    # parser.add_argument('-3b', '--3B', default='multijet/picoAOD_3b_wJCM_v0_newSBDef_2017_multijet.root')
    parser.add_argument('-s', '--seed', default=2435)
    parser.add_argument('-p', '--plot', default=False)
    args = parser.parse_args()
    main(args.seed, args.plot)

# import os
# os.system('python3 fakeData.py -s 7 -p True')#

# m4jPlot(m4jBinEdges, 
#             data3b = data3b, 
#             bgW = w3to4, #np.ones(len(data3b['m4j'])), 
#             dataMixed = data4b, 
#             dataW = np.ones(len(data4b['m4j'])), 
#             figName="test.png", plotAll=False, title='seed '+str(seed), m4jMax=1800, dijetLim = None, datLabel = None, show=True)
