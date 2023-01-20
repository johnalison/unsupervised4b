import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getDtoM_chunk, getDtoM                                                                                                                                                                     

parser = argparse.ArgumentParser()
# parser.add_argument('-d', '--data', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root")
# parser.add_argument('-t', '--ttbar', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root")
# parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_wkdt.root")
parser.add_argument('-m', '--mixed', default="1")
parser.add_argument('-s', '--selfCount', default="0")
parser.add_argument('-dYr', '--dataYear', default="2017")
parser.add_argument('-ttYr', '--ttYear', default="2017")
parser.add_argument('-ttP', '--ttProcess', default="TTTo2L2Nu")
parser.add_argument('-v', '--version', default="0")
parser.add_argument('-bf', '--binEdgeFile', default="unsupervised4b/m4jBinEdges.txt")#"../m4jBinEdges.txt")
# parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()


# print("Reading in data:\n",args.data)
# print("Reading in ttbar:\n",args.ttbar)

if args.mixed == "1":
    outputfile = "picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root"
    TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+args.ttProcess+args.ttYear+"_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root"
    dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+args.dataYear+"_3bDvTMix4bDvT_v"+args.version+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef.root"
else:
    outputfile = "picoAOD_3b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root" 
    TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+args.ttProcess+args.ttYear+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
    dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+args.dataYear+"_3b/picoAOD_3b_wJCM_newSBDef.root"

m4jBinEdges = np.array([[98,353],[353,406],[406,444],[444,480],[480,516],[516,558],[558,608],[608,679],[679,796],[796,1199]])
m4jBinEdges = np.loadtxt(args.binEdgeFile)
print(m4jBinEdges)
mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + args.version

def dataSel(data, binEdge):
    sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True,False))
    sel &= np.asarray(np.where((data['leadStM']>0) & (data['leadStM']<=250), True,False))
    sel &= np.asarray(np.where((data['sublStM']>0) & (data['sublStM']<=250), True,False))
    sel &= np.asarray(np.where(data['passHLT'], True, False))
    return sel


#
#  Reading in the branches
#

dataFile   = uproot.open(dataFilename)['Events']
arrayNames = [mcPsedoTagWeight,"nSelJets","event","run", "leadStM","sublStM","m4j","passHLT"]
arrayNames = [mcPsedoTagWeight,"passHLT","leadStM","sublStM","m4j"]
data       = dataFile.arrays(arrayNames)
print(arrayNames)

if args.selfCount == "0":
    ttbarFile = uproot.open(TTFilename)['Events']
    arrayNamesTTbar = [mcPsedoTagWeight,"nSelJets","event","run", "leadStM","sublStM","m4j","trigWeight_Data","passHLT"]
    arrayNamesTTbar = [mcPsedoTagWeight,"leadStM","sublStM","m4j","passHLT"]
    print(arrayNamesTTbar)
    ttbar      = ttbarFile.arrays(arrayNamesTTbar)
else:
    outputfile = (dataFilename.split("/")[-1]).split(".root")[0]+"_"+args.dataYear+"_wkdt.root"

print("Branches extracted")

wkdt = np.zeros(len(data['m4j']))
for binNo, binEdge in enumerate(m4jBinEdges):
    m4jPass = dataSel(data, binEdge)
    leadStM_data = data['leadStM'][m4jPass]
    sublStM_data = data['sublStM'][m4jPass]
    mcPsdoTagW = data[mcPsedoTagWeight][m4jPass]

    if args.selfCount == "0":
        m4jPassTT = dataSel(ttbar, binEdge)
        leadStM_ttbar = ttbar['leadStM'][m4jPassTT]
        sublStM_ttbar = ttbar['sublStM'][m4jPassTT]
        mcPsdoTagW_ttbar = ttbar[mcPsedoTagWeight][m4jPassTT]
        wkdt[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, radius = 5, ttWeight = mcPsdoTagW_ttbar)
    else:
        wkdt[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_data, sublStM_data, radius = 5, ttWeight = mcPsdoTagW)
    print("Bin selections obtained for bin", binNo)

print("DtoM reweighting complete")
with uproot.recreate(outputfile) as outfile:
    outfile["Events"] = {"wkdt": np.array(wkdt)}
    outfile["Events"].show()

print(outputfile, 'outputfile saved\n')