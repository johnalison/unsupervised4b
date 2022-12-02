import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getDtoM, getDtoM_chunk                                                                                                                                                                         

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root")
parser.add_argument('-t', '--ttbar', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root")
parser.add_argument('-pY', '--processYr', default="TTTo2L2Nu2017")
parser.add_argument('-v', '--version', default="0")
parser.add_argument('-bf', '--binEdgeFile', default="unsupervised4b/m4jBinEdges.txt")
parser.add_argument('-bin', '--bin', default=1)

parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_wDtoM.root")

#parser.add_argument('--kdTreeData', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX")
#parser.add_argument('--kdTreesTTbar', default="[root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX]")
parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()


print("Reading in data:\n",args.data)
print("Reading in ttbar:\n",args.ttbar)

dataFile   = uproot.open(args.data)['Events']
# ttbarFile = uproot.open(args.ttbar)['Events']
m4jBinEdges = np.loadtxt(args.binEdgeFile)

# print(ttbarFile.keys())


#
#  Reading in the branches
#

mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + args.version

arrayNames = [mcPsedoTagWeight,"nSelJets","event","run", "leadStM","sublStM","m4j","passHLT"]
arrayNames = [mcPsedoTagWeight,"passHLT","leadStM","sublStM","m4j"]
data       = dataFile.arrays(arrayNames)
print(arrayNames)
# data = {key.decode('utf-8'): value for (key, value) in data.items()}

# arrayNamesTTbar = [mcPsedoTagWeight,"nSelJets","event","run", "leadStM","sublStM","m4j","trigWeight_Data","passHLT"]
# arrayNamesTTbar = [mcPsedoTagWeight,"leadStM","sublStM","m4j","passHLT"]
# print(arrayNamesTTbar)
# ttbar      = ttbarFile.arrays(arrayNamesTTbar)

# ttbar = {key.decode('utf-8'): value for (key, value) in ttbar.items()}

print("Branches extracted")

def dataSel(data, binEdge):
    sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True,False))
    sel &= np.asarray(np.where((data['leadStM']>0) & (data['leadStM']<=250), True,False))
    sel &= np.asarray(np.where((data['sublStM']>0) & (data['sublStM']<=250), True,False))
    sel &= np.asarray(np.where(data['passHLT'], True, False))
    return sel

wDtoM = np.zeros(len(data['m4j']))
for binNo, binEdge in enumerate(m4jBinEdges):
    m4jPass = dataSel(data, binEdge)
    # m4jPassTT = dataSel(ttbar, binEdge)
    leadStM_data = data['leadStM'][m4jPass]
    sublStM_data = data['sublStM'][m4jPass]
    mcPsdoTagW = data[mcPsedoTagWeight][m4jPass]
    # leadStM_ttbar = ttbar['leadStM'][m4jPassTT]
    # sublStM_ttbar = ttbar['sublStM'][m4jPassTT]
    # mcPsdoTagW_ttbar = ttbar[mcPsedoTagWeight][m4jPassTT]
    # print("Bin selections obtained for bin", binNo)

    wDtoM[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_data, sublStM_data, radius = 5, ttWeight = mcPsdoTagW)
    # wDtoM[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, radius = 5, ttWeight = mcPsdoTagW_ttbar)

print("DtoM reweighting complete")
with uproot.recreate(args.outputfile) as outfile:
    outfile["Events"] = {"wDtoM": np.array(wDtoM)}
    outfile["Events"].show()


print(args.outputfile, 'outputfile saved\n')
