import argparse
import numpy as np
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda   
from kdtree import getDtoM, getDtoM_chunk                   
import warnings                                                                                                                                                      
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser()
# parser.add_argument('-4B', '--4b', default='multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0+_newSBDef_2017_multijet.root')
# parser.add_argument('-3b', '--3B', default='multijet/picoAOD_3b_wJCM_v0_newSBDef_2017_multijet.root')
parser.add_argument('-v', '--version', default=0)
parser.add_argument('-y', '--year', default="0")
parser.add_argument('-bf', '--binEdgeFile', default="../m4jBinEdges.txt")
# parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2018_w3to4.root")
# parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()

vn = args.version
year = args.year
m4jBinEdges = np.loadtxt(args.binEdgeFile)
# m4jBinEdges = np.loadtxt("../m4jBinEdges.txt")


outputfile = "w3to4_4b/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef_"+year+"_w3to4_4b.root"
outputfile = 'w3to4/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_w3to4.root'
wFilename4b = 'multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
filename4b = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+str(vn)+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef.root"

wFilename3b = 'multijet/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
filename3b = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"

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

arrayNames = ["passHLT","leadStM","sublStM","m4j"]
# print(arrayNames)

data4b = file4b.arrays(arrayNames)
data3b = file3b.arrays(arrayNames)
wDtoM4b = np.array(wFile4b.arrays('wDtoM'))['wDtoM']
wDtoM3b = np.array(wFile3b.arrays('wDtoM'))['wDtoM']
w3to4 = np.zeros(len(data3b['m4j']))
print(len(data4b['m4j']), len(data3b['m4j']), len(wDtoM4b), len(wDtoM3b))

def dataSel(data, binEdge, binEdgeSec = False):
    sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True, False))
    if binEdgeSec is not False:
        sel |= np.asarray(np.where((data['m4j']>binEdgeSec[0]) & (data['m4j']<=binEdgeSec[1]), True, False))
    sel &= np.asarray(np.where((data['leadStM']>0) & (data['leadStM']<=250), True,False))
    sel &= np.asarray(np.where((data['sublStM']>0) & (data['sublStM']<=250), True,False))
    sel &= np.asarray(np.where(data['passHLT'], True, False))
    return sel


for binNo in range(len(m4jBinEdges)):
    if binNo == 0 or binNo == len(m4jBinEdges)-1:
        continue
    m4jPass_SR3b = dataSel(data3b, m4jBinEdges[binNo])
    leadStM_SR3b = data3b['leadStM'][m4jPass_SR3b]
    sublStM_SR3b = data3b['sublStM'][m4jPass_SR3b]
    wDtoM_SR3b = wDtoM3b[m4jPass_SR3b]

    m4jPass_SB3b = dataSel(data3b, m4jBinEdges[binNo-1], m4jBinEdges[binNo+1])
    leadStM_SB3b = data3b['leadStM'][m4jPass_SB3b]
    sublStM_SB3b = data3b['sublStM'][m4jPass_SB3b]
    wDtoM_SB3b = wDtoM3b[m4jPass_SB3b]

    m4jPass_SB4b = dataSel(data4b, m4jBinEdges[binNo-1], m4jBinEdges[binNo+1])
    leadStM_SB4b = data4b['leadStM'][m4jPass_SB4b]
    sublStM_SB4b = data4b['sublStM'][m4jPass_SB4b]
    wDtoM_SB4b = wDtoM4b[m4jPass_SB4b]

    num = getDtoM(leadStM_SR3b, sublStM_SR3b, leadStM_SB4b, sublStM_SB4b, radius = 5, ttWeight = wDtoM_SB4b)
    den = getDtoM(leadStM_SR3b, sublStM_SR3b, leadStM_SB3b, sublStM_SB3b, radius = 5, ttWeight = wDtoM_SB3b)
    w3to4[m4jPass_SR3b] = np.divide(num,den)*wDtoM_SR3b

print("DtoM reweighting complete")
# np.savetxt(year+"_w3to4cor.txt", w3to4)
with uproot.recreate(outputfile) as outfile:
    outfile["Events"] = {"w3to4": np.array(w3to4)}
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
