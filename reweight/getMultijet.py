import numpy as np
import uproot

yearArr = ['2016', '2017', '2018']

processArr = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
#### mixed ####
for vn in range(0,15):
    for year in ['2017', '2018']:
        outputfile = 'multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
        dataFile = 'wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_wDtoM.root'
        mcpFile = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+str(vn)+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef.root"
        data = np.array(uproot.open(dataFile)['Events'].arrays('wDtoM'))['wDtoM']
        mcp = np.array(uproot.open(mcpFile)['Events'].arrays("mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)))["mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)]
        ttbar = 0
        for process in processArr:
            ttFile = 'wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+process+year+'_wDtoM.root'
            ttbar += np.array(uproot.open(ttFile)['Events'].arrays('wDtoM'))['wDtoM']
        with uproot.recreate(outputfile) as outfile:
            outfile["Events"] = {"wDtoM": np.divide(np.multiply(mcp,(data-ttbar)),data)}
            # outfile["Events"].show()
        print("Done", year, vn)

TTprocess = ['_preVFP','_postVFP']
for vn in range(15):
    for year in ['2016']:
        outputfile = 'multijet/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
        dataFile = 'wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_wDtoM.root'
        mcpFile = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+str(vn)+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+str(vn)+"_newSBDef.root"
        mcp = np.array(uproot.open(mcpFile)['Events'].arrays("mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)))["mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)]
        data = np.array(uproot.open(dataFile)['Events'].arrays('wDtoM'))['wDtoM']
        ttbar = 0
        for process in processArr:
            for tt in TTprocess:
                ttFile = 'wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v'+str(vn)+'_newSBDef_'+process+year+tt+'_wDtoM.root'
                ttbar += np.array(uproot.open(ttFile)['Events'].arrays('wDtoM'))['wDtoM']
        with uproot.recreate(outputfile) as outfile:
            outfile["Events"] = {"wDtoM": np.divide(np.multiply(mcp,(data-ttbar)),data)}
            # outfile["Events"].show()
        print("Done", year, vn)


#### 3b Data #####
import numpy as np
import uproot

yearArr = ['2016', '2017', '2018']

processArr = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']

for year in ['2017', '2018']:
    mcpFile = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
    mcpBranch = uproot.open(mcpFile)['Events']
    for vn in range(0,15):
        outputfile = 'multijet/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
        dataFile = "wDtoMSelf/picoAOD_3b_wJCM_v"+str(vn)+"_newSBDef_"+year+"_wDtoM.root"   
        data = np.array(uproot.open(dataFile)['Events'].arrays('wDtoM'))['wDtoM']
        mcp = np.array(mcpBranch.arrays("mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)))["mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)]
        ttbar = 0
        for process in processArr:
            ttFile = 'wDtoM3b/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+process+year+'_wDtoM.root'
            ttbar += np.array(uproot.open(ttFile)['Events'].arrays('wDtoM'))['wDtoM']
        with uproot.recreate(outputfile) as outfile:
            outfile["Events"] = {"wDtoM": np.divide(np.multiply(mcp,(data-ttbar)),data)}
            # outfile["Events"].show()
        print("Done", year, vn)

TTprocess = ['_preVFP','_postVFP']  
for year in ['2016']:
    mcpFile = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
    mcpBranch = uproot.open(mcpFile)['Events']
    for vn in range(15):
        outputfile = 'multijet/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+year+'_multijet.root'
        dataFile = "wDtoMSelf/picoAOD_3b_wJCM_v"+str(vn)+"_newSBDef_"+year+"_wDtoM.root"
        data = np.array(uproot.open(dataFile)['Events'].arrays('wDtoM'))['wDtoM']
        mcp = np.array(mcpBranch.arrays("mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)))["mcPseudoTagWeight_3bDvTMix4bDvT_v"+str(vn)]
        ttbar = 0
        for process in processArr:
            for tt in TTprocess:
                ttFile = 'wDtoM3b/picoAOD_3b_wJCM_v'+str(vn)+'_newSBDef_'+process+year+tt+'_wDtoM.root'
                ttbar += np.array(uproot.open(ttFile)['Events'].arrays('wDtoM'))['wDtoM']
        with uproot.recreate(outputfile) as outfile:
            outfile["Events"] = {"wDtoM": np.divide(np.multiply(mcp,(data-ttbar)),data)}
            # outfile["Events"].show()
        print("Done", year, vn)


### For Data ###
'''
wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2016_wDtoM.root
wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2017_wDtoM.root
wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_2018_wDtoM.root

wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTTo2L2Nu2016_postVFP_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTTo2L2Nu2016_preVFP_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToHadronic2016_postVFP_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToHadronic2016_preVFP_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToSemiLeptonic2016_postVFP_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToSemiLeptonic2016_preVFP_wDtoM.root

wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTTo2L2Nu2017_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToHadronic2017_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToSemiLeptonic2017_wDtoM.root

wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTTo2L2Nu2018_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToHadronic2018_wDtoM.root
wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_TTToSemiLeptonic2018_wDtoM.root
'''