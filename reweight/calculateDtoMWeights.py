import argparse
import uproot # https://github.com/scikit-hep/uproot3 is in lcg_99cuda                                                                                                                                                                      

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/mixed2016_3bDvTMix4bDvT_v0/picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef.root")
parser.add_argument('-o', '--outputfile', default="picoAOD_3bDvTMix4bDvT_4b_wJCM_v0_newSBDef_wDtoM.root")
parser.add_argument('-t', '--ttbar', default="[root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/ULTrig/TTTo2L2Nu2017_4b_noPSData_wTrigW/picoAOD_4b_wJCM_newSBDef.root,"])

#parser.add_argument('--kdTreeData', default="root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX")
#parser.add_argument('--kdTreesTTbar', default="[root://cmsxrootd.fnal.gov//store/user/jda102/condor/ZH4b/XXX]")
parser.add_argument('--debug', action="store_true", help='')
args = parser.parse_args()


print("Reading in data",args.data)
#print("Reading in ttar",args.ttbar)
dataFile=uproot.open(args.data)['Events']
#ttbarFile=uproot.open(args.ttbar)['Events']

#
#  Reading in the branches
#
arrayNames = ["mcPseudoTagWeight_3bDvTMix4bDvT_v0","nSelJets","event","run", "leadStM","m4j","sublStM","passHLT","trigWeight_Data"]
data       = dataFile.arrays(arrayNames)

arrayNamesTTbar = ["mcPseudoTagWeight_3bDvTMix4bDvT_v0"]
ttbar      = ttbarFile.arrays(arrayNamesTTbar)



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

    outfile["Events"] = {"wDToM": np.array(wDtM), "event": np.array(data["event"][:maxEvents]), "run":np.array(data["run"][:maxEvents])}
    outfile["Events"].show()

writeout (wDtoM)


#dataM4jSR     = function(m4jlow, m4jhigh)
#dataM4jSB_up  = function(m4jlow, m4jhigh)
#dataM4jSB_low = function(m4jlow, m4jhigh)
#
#doAnalysis(dataM4jS, dataM4jSB_up, dataM4jSB_low)
