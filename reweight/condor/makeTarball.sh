localPathToTarball=/uscms/home/jda102/nobackup/HH4b/semisupervised
CMSSW=CMSSW_11_1_0_pre5
rm -rf ${localPathToTarball}/${CMSSW}.tgz
cd ${localPathToTarball}
tar -C ${localPathToTarball} -zcvf ${localPathToTarball}/${CMSSW}.tgz ${CMSSW} --exclude="*root" --exclude-vcs --exclude-caches-all
ls ${localPathToTarball} -alh
xrdfs root://cmseos.fnal.gov/ mkdir /store/user/jda102/condor
xrdcp -f ${localPathToTarball}/${CMSSW}.tgz root://cmseos.fnal.gov//store/user/jda102/condor/${CMSSW}.tgz
cd -
