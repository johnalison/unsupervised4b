# from fileinput import filename
import numpy as np; import uproot; import matplotlib.pyplot as plt;

from data_class import file

print("\n")

# filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2016_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root"
# dataFile = file(name = filename, type = "mix")
# dataFile.getSelections()
# print("done", "\n")
# del dataFile

# filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data2016_3b/picoAOD_3b_wJCM_newSBDef.root"
# dataFile = file(name = filename, type = "3b")
# dataFile.getSelections()
# print("done", "\n")
# del dataFile

#filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
##filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root"
#ttFile = file(name = filename, type = "3b")
#ttFile.getTTSelections()
#print("done", "\n")
#del ttFile
#

val = []
processArr = ["TTTo2L2Nu", "TTToHadronic", "TTToSemiLeptonic"]
# yearArr = ["2016_postVFP","2016_preVFP"]
yearArr = ['2018']
for year in yearArr:
    for process in processArr:
        #filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+process+year+"_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root"
        #ttFile = file(name = filename, type = "4btt")
        filename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+process+year+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
        ttFile = file(name = filename, type = "4btt")
        # fttFile.getTTSelections()
        val.append(ttFile.selTTAll(True))

allprocess = np.sum(val, axis = 0)
print("\n", "Counts:", int(allprocess[0]), "    Weights:", allprocess[1] )
print("done", "\n")
