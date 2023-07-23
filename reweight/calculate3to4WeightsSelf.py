import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getkNNwJCM, getDtoMwJCM               
from kdtree import dataSel
import awkward as ak
import warnings                                                                                                                                                      
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser()
# parser.add_argument('-4B', '--4b', default='multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0+_newSBDef_2017_multijet.root')
# parser.add_argument('-3b', '--3B', default='multijet/picoAOD_3b_wJCM_v0_newSBDef_2017_multijet.root')
parser.add_argument('-v', '--version', default=0)
parser.add_argument('-w', '--weightType', default='weight')
parser.add_argument('-y', '--year', default="2017")
parser.add_argument('-sig', '--sigIn', default="0")
parser.add_argument('-bf', '--binEdgeFile', default="unsupervised4b/m4jBinEdges.txt")
parser.add_argument('-mean', '--mean', default='800')
parser.add_argument('-std', '--std', default='30')
parser.add_argument('-p', '--pathDtoM', default="root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/wDtoMwJCM/")
# parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2018_w3to4.root")
# parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()

vn = args.version
year = args.year
# m4jBinEdges = np.loadtxt(args.binEdgeFile)
# m4jBinEdges = np.loadtxt("../m4jBinEdges.txt")
# m4jBinEdges = np.array([[98,353],[353,406],[406,444],[444,480],[480,516],[516,558],[558,608],[608,679],[679,796],[796,1199]])
m4jBinEdges = np.array([[0, 361], [361, 425], [425, 479], [479, 533], [533, 591], [591, 658], [658, 741], [741, 854], [854, 1044], [1044, 1800]])

def gaus(x, mu, sigma):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))


# outputfile = 'w3to4/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_w3to4.root'
# wFilename4b = 'multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
# filename4b = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+str(vn)+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef.root"
# wFilename3b = 'multijet/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
# filename3b = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"

EOSwDtoMpath = args.pathDtoM

wFilename4b = EOSwDtoMpath+"mixed"+year+"_picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef_wDtoM.root"
filename4b = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/mixed"+year+"_picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef.root"
wFilename3b = EOSwDtoMpath+"data"+year+"_picoAOD_3b_wJCM_newSBDef_wDtoM.root"
filename3b = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"data"+year+"_picoAOD_3b_wJCM_newSBDef.root"
outputfileL =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_w3to4SBL.root"
outputfileU =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_w3to4SBU.root"
outputfile =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_w3to4SR.root"


sigFilename = None
if args.sigIn != '0':
    wFilename4b = EOSwDtoMpath+"mixed"+year+"_picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef_ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_wDtoM.root"
    outputfileL =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_w3to4SBL.root"
    outputfileU =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_w3to4SBU.root"
    outputfile =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_ZH4b_"+str(args.sigIn)+"percentGauss"+str(args.mean)+"-"+str(args.std)+"_w3to4SR.root"
    sigFilename = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"ZH4b"+year+"_picoAOD.root"

# filename4b = "root://cmseos.fnal.gov//store/user/smurthy/condor/unsupervised4b/randPair/files/"+"ZH4b"+year+"_picoAOD.root"
# outputfile =  "data"+year+"_picoAOD_3b_wJCM_newSBDef_ZH4b_w3to4.root"

print(vn, year)
print("Reading in files:")
print(wFilename4b)
print(filename4b)
print(wFilename3b)
print(filename3b)


file4b = uproot.open(filename4b)['Events']
wFile4b = uproot.open(wFilename4b)['Events']

file3b = uproot.open(filename3b)['Events']
wFile3b = uproot.open(wFilename3b)['Events']
mcPseudoTagWeight = "mcPseudoTagWeight"#_3bDvTMix4bDvT_v" + args.version
mcPseudoTagWeight = 'weight'
mcPseudoTagWeight = args.weightType

arrayNames = ["passHLT","leadStM","sublStM","m4j",mcPseudoTagWeight,'nSelJets','SR']
# print(arrayNames)

data4b = file4b.arrays(arrayNames)

if sigFilename is not None:
    sigFile   = uproot.open(sigFilename)['Events']
    sig       = sigFile.arrays(arrayNames)
    m4jPass = dataSel(data4b,[m4jBinEdges[0][0], m4jBinEdges[-1][1]])
    mcPsdoTagWSum = np.sum(np.array(data4b[mcPseudoTagWeight][m4jPass]),dtype=np.float32)
    sigPass = dataSel(sig, [m4jBinEdges[0][0], m4jBinEdges[-1][1]])
    sigWSum = np.sum(np.array(sig[mcPseudoTagWeight][sigPass]),dtype=np.float32)
    sigFactor = float(args.sigIn)*0.01*mcPsdoTagWSum
    sigRescale = np.zeros(len(sigPass))
    gausScale = gaus(sig['m4j'][sigPass], mu=float(args.mean), sigma=float(args.std))
    sigRescale[sigPass] = sigFactor * gausScale/np.sum(np.array(gausScale), dtype=np.float32)
    # sigRescale[sigPass] = sigFactor * sig[mcPseudoTagWeight][sigPass]/sigWSum
    sig[mcPseudoTagWeight] = sigRescale
    data4b = ak.concatenate((data4b, sig))

data3b = file3b.arrays(arrayNames)
wDtoM4b = np.multiply(np.array(wFile4b.arrays('wDtoM'))['wDtoM'],data4b[mcPseudoTagWeight])
wDtoM3b = np.multiply(np.array(wFile3b.arrays('wDtoM'))['wDtoM'],data3b[mcPseudoTagWeight])

w3to4_L = np.zeros(len(data3b['m4j']))
w3to4_U = np.zeros(len(data3b['m4j']))
w3to4SR = np.zeros(len(data3b['m4j']))
print(len(data4b['m4j']), len(data3b['m4j']), len(wDtoM4b), len(wDtoM3b))

for binNo in range(len(m4jBinEdges)):
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
    nSelJets_SB3b = data3b['nSelJets'][m4jPass_SB3b]
    wDtoM_SB3b = wDtoM3b[m4jPass_SB3b]
    binNoSB3b = np.array(data3b['SR'][m4jPass_SB3b])
    # wDtoM_SB3b[np.where(np.isnan(wDtoM_SB3b)==True)] = 0

    leadStM_SB4b = data4b['leadStM'][m4jPass_SB4b]
    sublStM_SB4b = data4b['sublStM'][m4jPass_SB4b]
    nSelJets_SB4b = data4b['nSelJets'][m4jPass_SB4b]
    wDtoM_SB4b = wDtoM4b[m4jPass_SB4b]
    # wDtoM_SB4b[np.where(np.isnan(wDtoM_SB4b)==True)] = 0

    m4jPass_SR3b = dataSel(data3b, m4jBinEdges[binNo])
    leadStM_SR3b = data3b['leadStM'][m4jPass_SR3b]
    sublStM_SR3b = data3b['sublStM'][m4jPass_SR3b]
    nSelJets_SR3b = data3b['nSelJets'][m4jPass_SR3b]
    wDtoM_SR3b = wDtoM3b[m4jPass_SR3b]
    # wDtoM_SR3b[np.where(np.isnan(wDtoM_SR3b)==True)] = 0

    w3to4 = np.zeros(len(leadStM_SB3b))
    NN = 35
    if mcPseudoTagWeight == "weight":
        # num = getDtoMwJCM(leadStM_SB3b, sublStM_SB3b, leadStM_SB4b, sublStM_SB4b, radius = 15, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SB3b, nSJtt = nSelJets_SB4b)
        # den = getDtoMwJCM(leadStM_SB3b, sublStM_SB3b, leadStM_SB3b, sublStM_SB3b, radius = 15, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SB3b, nSJtt = nSelJets_SB3b)

        # numSR = getDtoMwJCM(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, radius = 15, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB4b)
        # denSR = getDtoMwJCM(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, radius = 15, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB3b)


        num = getkNNwJCM(leadStM_SB3b, sublStM_SB3b, leadStM_SB4b, sublStM_SB4b, NN=NN, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SB3b, nSJtt = nSelJets_SB4b)
        den = getkNNwJCM(leadStM_SB3b, sublStM_SB3b, leadStM_SB3b, sublStM_SB3b, NN=NN, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SB3b, nSJtt = nSelJets_SB3b)
        
        numSR = getkNNwJCM(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, NN=NN, ttWeight = wDtoM_SB4b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB4b)
        denSR = getkNNwJCM(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, NN=NN, ttWeight = wDtoM_SB3b, nSJDat = nSelJets_SR3b, nSJtt = nSelJets_SB3b)

    w3to4SR[m4jPass_SR3b] = np.divide(numSR,denSR)

    w3to4 = np.divide(num,den)
    if binNo == 0:
        w3to4_U[dataSel(data3b, m4jBinEdges[binNo+1])] = w3to4[binNoSB3b==binNo+1]
    elif binNo == len(m4jBinEdges)-1:
        w3to4_L[dataSel(data3b, m4jBinEdges[binNo-1])] = w3to4[binNoSB3b==binNo-1]
    else:
        w3to4_L[dataSel(data3b, m4jBinEdges[binNo-1])] = w3to4[binNoSB3b==binNo-1]
        w3to4_U[dataSel(data3b, m4jBinEdges[binNo+1])] = w3to4[binNoSB3b==binNo+1]
    

print(w3to4_U)
print(w3to4_L)
print(w3to4SR)
print("3to4 reweighting complete")
# np.savetxt(year+"_w3to4cor.txt", w3to4)

with uproot.recreate(outputfileL) as outfile:
    outfile["Events"] = {"w3to4": np.array(w3to4_L)}
    # outfile["Events"].show()
print(outputfileL, 'outputfile saved\n')

with uproot.recreate(outputfileU) as outfile:
    outfile["Events"] = {"w3to4": np.array(w3to4_U)}
    # outfile["Events"].show()
print(outputfileU, 'outputfile saved\n')

with uproot.recreate(outputfile) as outfile:
    outfile["Events"] = {"w3to4": np.array(w3to4SR)}
    # outfile["Events"].show()
print(outputfile, 'outputfile saved\n')


'''
all for one year 2018
get 3b and 4b
get 3b and 4b wDtoM 
bin data based on binbounds
w3to4 = 3bSRwDtoM*kdtree(3bSR, 4bSB, wDtoM)/kdtree(3bSR, 3bSB, wDtoM)
will have same dimensions as 3b data
save file
'''
