import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getKDEwJMC #getkNNwJCM#getDtoM_chunk, getDtoMwJCM     
from kdtree import dataSel
import awkward as ak                                                                                                                                                         

parser = argparse.ArgumentParser()
# parser.add_argument('-d', '--data', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root")
# parser.add_argument('-t', '--ttbar', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root")
# parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_wkdt.root")
parser.add_argument('-m', '--mixed', default="m")
parser.add_argument('-s', '--selfCount', default="0")
parser.add_argument('-dYr', '--dataYear', default="2016")
parser.add_argument('-ttYr', '--ttYear', default="2016_postVFP")
parser.add_argument('-ttP', '--ttProcess', default="TTToHadronic")
parser.add_argument('-v', '--version', default="0")
parser.add_argument('-sig', '--sigIn', default="0")
parser.add_argument('-w', '--weightType', default='weight')
parser.add_argument('-mean', '--mean', default='800')
parser.add_argument('-std', '--std', default='30')
parser.add_argument('-bf', '--binEdgeFile', default="unsupervised4b/m4jBinEdges.txt")#"../m4jBinEdges.txt")
# parser.add_argument('--debug', action="store_true", help='')

args = parser.parse_args()


# print("Reading in data:\n",args.data)
# print("Reading in ttbar:\n",args.ttbar)

# if args.mixed == "m":
#     outputfile = "picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root"
#     TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+args.ttProcess+args.ttYear+"_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root"
#     dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+args.dataYear+"_3bDvTMix4bDvT_v"+args.version+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef.root"
# else:
#     outputfile = "picoAOD_3b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root" 
#     TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+args.ttProcess+args.ttYear+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#     dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+args.dataYear+"_3b/picoAOD_3b_wJCM_newSBDef.root"

sigFilename = None
if args.mixed == "m":
    outputfile = "picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root"
    TTFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+args.ttProcess+args.ttYear+"_4b_noPSData_wTrigW_picoAOD_4b_wJCM_newSBDef.root"
    dataFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"mixed"+args.dataYear+"_picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+args.version+"_newSBDef.root"
    if args.sigIn != '0':
        outputfile = outputfile.split("wkdt.root")[0]+"ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_wkdt.root"
        sigFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"ZH4b"+args.dataYear+"_picoAOD.root"
elif args.mixed == 'd':
    outputfile = "picoAOD_3b_wJCM_v"+args.version+"_newSBDef_"+args.ttProcess+args.ttYear+"_wkdt.root" 
    TTFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+args.ttProcess+args.ttYear+"_3b_wTrigW_picoAOD_3b_wJCM_newSBDef.root"
    dataFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"data"+args.dataYear+"_picoAOD_3b_wJCM_newSBDef.root"

print(outputfile)
print(TTFilename)
print(dataFilename)
print(sigFilename)

mixEventCount = [115186.0, 111388.0, 161797.0]
def gaus(x, mu, sigma):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))

# m4jBinEdges = np.loadtxt(args.binEdgeFile)
# m4jBinEdges = np.array([[98,353],[353,406],[406,444],[444,480],[480,516],[516,558],[558,608],[608,679],[679,796],[796,1199]])
m4jBinEdges = np.array([[0, 361], [361, 425], [425, 479], [479, 533], [533, 591], [591, 658], [658, 741], [741, 854], [854, 1044], [1044, 1800]])

print(m4jBinEdges)
# mcPseudoTagWeight = "mcPseudoTagWeightInclusive"#_3bDvTMix4bDvT_v" + args.version
# mcPseudoTagWeight = 'weight'
# mcPseudoTagWeight = "mcPseudoTagWeight"
mcPseudoTagWeight = args.weightType
arrayNames = [mcPseudoTagWeight,"passHLT","leadStM","sublStM","m4j", 'nSelJets']
print(arrayNames)
#
#  Reading in the branches
#

dataFile   = uproot.open(dataFilename)['Events']
data       = dataFile.arrays(arrayNames)

if sigFilename:
    sigFile   = uproot.open(sigFilename)['Events']
    sig       = sigFile.arrays(arrayNames)
    m4jPass = dataSel(data,[m4jBinEdges[0][0], m4jBinEdges[-1][1]])
    mcPsdoTagWSum = np.sum(np.array(data[mcPseudoTagWeight][m4jPass]),dtype=np.float32)
    sigPass = dataSel(sig, [m4jBinEdges[0][0], m4jBinEdges[-1][1]])
    sigWSum = np.sum(np.array(sig[mcPseudoTagWeight][sigPass]),dtype=np.float32)
    sigFactor = float(args.sigIn)*0.01*mcPsdoTagWSum
    gausScale = gaus(sig['m4j'][sigPass], mu=float(args.mean), sigma=float(args.std))
    sigRescale = np.zeros(len(sigPass))
    # sigRescale[sigPass] = sigFactor * sig[mcPseudoTagWeight][sigPass]/sigWSum
    sigRescale[sigPass] = sigFactor * gausScale/np.sum(np.array(gausScale), dtype=np.float32)
    sig[mcPseudoTagWeight] = sigRescale
    data = ak.concatenate((data, sig))

if args.selfCount == "0":
    ttbarFile = uproot.open(TTFilename)['Events']
    arrayNamesTTbar = arrayNames
    print(arrayNamesTTbar)
    ttbar      = ttbarFile.arrays(arrayNamesTTbar)
else:
    outputfile = (dataFilename.split("/")[-1]).split(".root")[0]+"_wkdt.root"
    if args.sigIn != '0' :
        outputfile = (dataFilename.split("/")[-1]).split(".root")[0]+"_ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_wkdt.root"

print("Branches extracted")

wkdt = np.zeros(len(data['m4j']))
for binNo, binEdge in enumerate(m4jBinEdges):
    m4jPass = dataSel(data, binEdge)
    leadStM_data = data['leadStM'][m4jPass]
    sublStM_data = data['sublStM'][m4jPass]
    nSelJets_data = data['nSelJets'][m4jPass]
    mcPsdoTagW = data[mcPseudoTagWeight][m4jPass]

    sigma = 5
    NN=35
    if args.selfCount == "0":
        m4jPassTT = dataSel(ttbar, binEdge)
        leadStM_ttbar = ttbar['leadStM'][m4jPassTT]
        sublStM_ttbar = ttbar['sublStM'][m4jPassTT]
        nSelJets_ttbar = ttbar['nSelJets'][m4jPassTT]
        mcPsdoTagW_ttbar = ttbar[mcPseudoTagWeight][m4jPassTT]
        wkdt[m4jPass] = getKDEwJMC(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, ttWeight = mcPsdoTagW_ttbar, nSJDat = nSelJets_data, nSJtt = nSelJets_ttbar, sigma = sigma)
        # wkdt[m4jPass] = getkNNwJCM(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, NN=NN, ttWeight = mcPsdoTagW_ttbar, nSJDat = nSelJets_data, nSJtt = nSelJets_ttbar)
        # wkdt[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, radius = 5, ttWeight = mcPsdoTagW_ttbar)
        # wkdt[m4jPass] = getDtoMwJCM(leadStM_data, sublStM_data, leadStM_ttbar, sublStM_ttbar, radius = 15, ttWeight = mcPsdoTagW_ttbar, nSJDat = nSelJets_data, nSJtt = nSelJets_ttbar)
    else:
        wkdt[m4jPass] = getKDEwJMC(leadStM_data, sublStM_data, leadStM_data, sublStM_data, ttWeight = mcPsdoTagW, nSJDat = nSelJets_data, nSJtt = nSelJets_data, sigma = sigma)
        # wkdt[m4jPass] = getDtoM_chunk(leadStM_data, sublStM_data, leadStM_data, sublStM_data, radius = 5, ttWeight = mcPsdoTagW)
        # wkdt[m4jPass] = getDtoMwJCM(leadStM_data, sublStM_data, leadStM_data, sublStM_data, radius = 15, ttWeight = mcPsdoTagW, nSJDat = nSelJets_data, nSJtt = nSelJets_data)
    print("Bin selections obtained for bin", binNo)

print("DtoM reweighting complete")
with uproot.recreate(outputfile) as outfile:
    outfile["Events"] = {"wkdt": np.array(wkdt)}
    outfile["Events"].show()

print(outputfile, 'outputfile saved\n')