import os 
import time

start = time.time()
calcDtoMPath = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights.py'
# calcDtoMPath3b = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights3b.py'

# processArr = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
# TTprocess = ['_preVFP','_postVFP','','']
# yearArr = ['2016','2016','2017','2018']
yearArr = ['2016','2017','2018']

vx = [str(i) for i in range(0,15)]

for vn in vx:
    for year in (yearArr):    
        outfile = "wDtoMSelf/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+vn+"_newSBDef_"+year+"_wDtoM.root"
        dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+vn+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+vn+"_newSBDef.root"
        os.system('python3 ' + calcDtoMPath + ' -v ' + vn + ' -d ' + dataFilename + ' -o ' +outfile) 


## Mixed Data ###
# for vn in vx:
#     for y, year in enumerate(yearArr):
#         for process in processArr:
#             processYear = process+year
#             print(processYear, vn)
#             outfile = "wDtoM/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+vn+"_newSBDef_"+processYear+TTprocess[y]+"_wDtoM.root"
#             TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+TTprocess[y]+"_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root"
#             dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed"+year+"_3bDvTMix4bDvT_v"+vn+"/picoAOD_3bDvTMix4bDvT_4b_wJCM_v"+vn+"_newSBDef.root"
#             os.system('python3 ' + calcDtoMPath + ' -pY ' + processYear + ' -v ' + vn + ' -d ' + dataFilename + ' -t ' + TTFilename + ' -o ' +outfile) 

# ### 3b Data ###

# for y, year in enumerate(yearArr):
#     for process in processArr:
#         processYear = process+year
#         print(processYear)
#         # outfile = str(["wDtoM3b/picoAOD_3b_wJCM_v"+vn+"_newSBDef_"+processYear+TTprocess[y]+"_wDtoM.root" for vn in range(0,15)])
#         TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+TTprocess[y]+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#         dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
#         os.system('python3 ' + calcDtoMPath3b + ' -pY ' + processYear+ ' -TTP ' + TTprocess[y] + ' -d ' + dataFilename + ' -t ' + TTFilename) 


# end = time.time()
# print("Done!", end - start)


# import os 
# import time
# start = time.time()
# calcDtoMPath = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights.py'
 
# processArr = ['TTToHadronic', 'TTToSemiLeptonic']
# TTprocess = ['_preVFP']
# yearArr = ['2016']
 
# vx = [str(i) for i in range(12,15)]
 
# for vn in vx:
#     for y, year in enumerate(yearArr):
#         for process in processArr:
#             processYear = process+year
#             print(processYear, vn)
#             outfile = "wDtoM3b/picoAOD_3b_wJCM_v"+vn+"_newSBDef_"+processYear+TTprocess[y]+"_wDtoM.root" 
#             TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+TTprocess[y]+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#             dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
#             os.system('python3 ' + calcDtoMPath + ' -pY ' + processYear+ ' -d ' + dataFilename + ' -t ' + TTFilename + ' -o ' + outfile + ' -v ' + vn)  



# import os 
# import time
# start = time.time()
# calcDtoMPath = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights.py'
 
# processArr = ['TTToHadronic', 'TTToSemiLeptonic']
# TTprocess = ['_preVFP']
# yearArr = ['2016']
 
# vx = [str(i) for i in range(12,15)]
 
# for vn in vx:
#     for y, year in enumerate(yearArr):
#         for process in processArr:
#             processYear = process+year
#             print(processYear, vn)
#             outfile = "wDtoM3b/picoAOD_3b_wJCM_v"+vn+"_newSBDef_"+processYear+TTprocess[y]+"_wDtoM.root" 
#             TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+TTprocess[y]+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#             dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
#             os.system('python3 ' + calcDtoMPath + ' -pY ' + processYear+ ' -d ' + dataFilename + ' -t ' + TTFilename + ' -o ' + outfile + ' -v ' + vn)  



# import os 
# import time
# start = time.time()
# calcDtoMPath = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights.py'
 
# processArr = ['TTToHadronic']
# # processArr = ['TTToSemiLeptonic']
# yearArr = ['2018']
 
# vx = ['1','3','5','7','10','12','13','14']

# vx = ['3']
 
# for vn in vx:
#     for y, year in enumerate(yearArr):
#         for process in processArr:
#             processYear = process+year
#             print(processYear, vn)
#             outfile = "wDtoM2/picoAOD_3b_wJCM_v"+vn+"_newSBDef_"+processYear+"_wDtoM.root" 
#             TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#             dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
#             os.system('python3 ' + calcDtoMPath + ' -pY ' + processYear+ ' -d ' + dataFilename + ' -t ' + TTFilename + ' -o ' + outfile + ' -v ' + vn)  

# print("done")

# import os 
# import time
# start = time.time()
# calcDtoMPath = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculateDtoMWeights.py'
 

# yearArr = ['2018']

# vx = [str(i) for i in range(12,15)]
 
# for vn in vx:
#     for y, year in enumerate(yearArr):
#         for process in processArr:
#             processYear = process+year
#             print(processYear, vn)
#             outfile = "wDtoM2/picoAOD_3b_wJCM_v"+vn+"_newSBDef_"+processYear+"_wDtoM.root" 
#             TTFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/"+processYear+"_3b_wTrigW/picoAOD_3b_wJCM_newSBDef.root"
#             dataFilename = "root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/data"+year+"_3b/picoAOD_3b_wJCM_newSBDef.root"
#             os.system('python3 ' + calcDtoMPath + ' -pY ' + processYear+ ' -d ' + dataFilename + ' -t ' + TTFilename + ' -o ' + outfile + ' -v ' + vn)  

# print("done")