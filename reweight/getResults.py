import numpy as np; import uproot
from getQuantilePlots import getQuantile
from getQuantilePlots import getDataForPlot, m4jPlot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

m4jBinEdges = np.loadtxt("../m4jBinEdges.txt")

# filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root"
# dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root"
# w3to4Filename = "w3to4_4b/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2017_w3to4_4b.root"

filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data2017_3b/picoAOD_3b_wJCM_newSBDef.root"
dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2017_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root"
w3to4Filename = "w3to4/picoAOD_3b_wJCM_v0_newSBDef_2017_w3to4.root"


processArr = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
#change year
vn = 0
ttFilename = ["root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+p+"2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root" for p in processArr]


####################################################

data3b, bgW, dataMixed = getDataForPlot(filename, dataFilename, w3to4Filename, ttFilename = None, vn = 0)
print(len(data3b['m4j']), len(bgW))
print("Data obtained")

m4jPlot(m4jBinEdges, data3b, bgW, dataMixed, figName = "m4j3to4DtoM.png")
exit()

####################################################

pp = PdfPages('foot.pdf')
# for binNo, m4jBinEdge in enumerate(m4jBinEdges[1:2]):
for binNo, m4jBinEdge in enumerate(m4jBinEdges):
    if binNo == 0 or binNo == len(m4jBinEdges)-1:
        continue

    m4jPass_SR3b = getQuantile.dataSel(data3b, m4jBinEdges[binNo])
    leadStM_SR3b = data3b['leadStM'][m4jPass_SR3b]
    sublStM_SR3b = data3b['sublStM'][m4jPass_SR3b]
    w3b = bgW[m4jPass_SR3b]
    
    m4jPass_SRMixed = getQuantile.dataSel(dataMixed, m4jBinEdges[binNo])
    leadStM_mixed = dataMixed['leadStM'][m4jPass_SRMixed]
    sublStM_mixed = dataMixed['sublStM'][m4jPass_SRMixed]
    wMixed = dataMixed['mcPseudoTagWeight_3bDvTMix4bDvT_v0'][m4jPass_SRMixed]

    qData = getQuantile(leadStM_SR3b, sublStM_SR3b, threshold = 35, weights = w3b, binNo=binNo, m4jBinEdge=m4jBinEdge)
    qData.defData(leadStM_mixed, sublStM_mixed, dataW = wMixed)
    figf = qData.getHistPlots()
    
    pp.savefig(figf); 
    plt.show();
    plt.close()

pp.close()

print("Done!")



 
