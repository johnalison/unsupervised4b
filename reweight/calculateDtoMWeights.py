import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getDtoM                                                                                                                                                                   

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root")
parser.add_argument('-t', '--ttbar', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root")

parser.add_argument('-bf', '--binEdgeFile', default="../m4jBinEdges.txt")
parser.add_argument('-bin', '--bin', default=1)

parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_wDtoM.root")

#parser.add_argument('--kdTreeData', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX")
#parser.add_argument('--kdTreesTTbar', default="[root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX]")
parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()


print("Reading in data:\n",args.data)
print("Reading in ttbar:\n",args.ttbar)

dataFile   = uproot.open(args.data)['Events']
ttbarFile = uproot.open(args.ttbar)['Events']
m4jBinEdges = np.loadtxt(args.binEdgeFile)


#
#  Reading in the branches
#

arrayNames = ["mcPseudoTagWeight_3bDvTMix4bDvT_v0","passHLT"]
arrayNames = ["mcPseudoTagWeight_3bDvTMix4bDvT_v0","nSelJets","event","run", "leadStM","sublStM","m4j","passHLT"]
data       = dataFile.arrays(arrayNames)
data = {key.decode('utf-8'): value for (key, value) in data.items()}

arrayNamesTTbar = ["mcPseudoTagWeight_3bDvTMix4bDvT_v0","nSelJets","event","run", "leadStM","sublStM","m4j","trigWeight_Data","passHLT"]
ttbar      = ttbarFile.arrays(arrayNamesTTbar)
ttbar = {key.decode('utf-8'): value for (key, value) in ttbar.items()}

print("Branches extracted")

def dataSel(data, binEdge):
    sel = np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True,False)
    sel &= np.where((data['leadStM']>0) & (data['leadStM']<=250), True,False)
    sel &= np.where((data['sublStM']>0) & (data['sublStM']<=250), True,False)
    sel &= np.where(data['passHLT'], True, False)
    return sel

wDtoM = np.zeros(len(data['m4j']))


# for binNo, binEdge in enumerate(m4jBinEdges):
binNo, binEdge = 0, m4jBinEdges[0]
# m4jPass = np.logical_and(data['m4j'] >= binEdge[0], data['m4j'] < binEdge[1])
m4jPass = dataSel(data, binEdge)
m4jPassTT = dataSel(ttbar, binEdge)
leadStM_data = data['leadStM'][m4jPass]
sublStM_data = data['sublStM'][m4jPass]
mcPsdoTagW = data['mcPseudoTagWeight_3bDvTMix4bDvT_v0'][m4jPass]
leadStM_ttbar = ttbar['leadStM'][m4jPassTT]
sublStM_ttbar = ttbar['sublStM'][m4jPassTT]
mcPsdoTagW_ttbar = ttbar['mcPseudoTagWeight_3bDvTMix4bDvT_v0'][m4jPassTT]
print("Bin selections obtained")

wDtoM[m4jPass] = getDtoM(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, radius = 5, ttWeight = mcPsdoTagW_ttbar)
print("DtoM reweighting complete")

# with uproot.recreate(args.outputfile) as outfile:
#     outfile = {"wDtoM": np.array(wDtoM)}
    # outfile["Events"].show()
with uproot.recreate(args.outputfile) as outfile:
    outfile["Events"] = {"wDtoM": np.array(wDtoM)}
    outfile["Events"].show()

print('outputfile saved')


# outfile = uproot.recreate(args.outputfile)
'''
loop on m4j bins 

    filter = data[:]["m4j"] > xx &&  data[:]["m4j"] < xx
    dataFile
    select data  in that bin     
    select ttbar  in that bin     

    Make kdTree data 
    Make kdTrees ttbar

    loop on data in that bin

        #
        # find indices of NN in data
        # 
        totalWeightData = 0
        arr_ind = getNN(data.leadStM, data.sublStM)
        loop on ind: totalWeightData += weight[ind]
            

        loop on ttbar kdTree
           find nNN in this ttbar kdTree
           See above

        compute weights
        addweight to weights_per_event_per_m4j_bin


#
#  Loop to convert wDToM[m4jBin] => wDToM[event]
#
indicesRead = [] 
wDtoM = []
loop on data
   find m4j bin
   look up weight in weights_per_event_per_m4j_bin
   wDtoM.append(weight)
   increment the indices read


with uproot.recreate(args.output) as outfile:

    outfile["Events"] = {"wDtoM": np.array(wDtoM), "event": np.array(data["event"][:maxEvents]), "run":np.array(data["run"][:maxEvents])}
    outfile["Events"].show()

writeout (wDtoM)


#dataM4jSR     = function(m4jlow, m4jhigh)
#dataM4jSB_up  = function(m4jlow, m4jhigh)
#dataM4jSB_low = function(m4jlow, m4jhigh)
#
#doAnalysis(dataM4jS, dataM4jSB_up, dataM4jSB_low)
'''