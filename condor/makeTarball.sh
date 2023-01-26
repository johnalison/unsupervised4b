CMSSWpath=/uscms/home/smurthy/nobackup/
localPathToTarball=/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b
CMSSW=CMSSW_11_1_0_pre5
# Remove all files in hierarchy and don't ask for confirmation
# The tgz file is located at that path with that name. Remove tarball
rm -rf ${localPathToTarball}/${CMSSW}.tgz
# zcvf = create archive and zip with given filename and display verbose; -C in specified directory
# exclude any root files and version-control-files and cache
# All files in CMSSW are archived into the .tgz file
cd ${localPathToTarball}
tar -C ${CMSSWpath} -zcvf ${localPathToTarball}/${CMSSW}.tgz ${CMSSW} --exclude="*root" --exclude="*npy" --exclude="*.pdf" --exclude="*.jdl" --exclude="*.stdout" --exclude="*.stderr" --exclude="*.log" --exclude="*.png" --exclude="*.dat" --exclude="*.tgz" --exclude-vcs --exclude-caches-all 
# All long-list human readable
ls ${localPathToTarball} -alh
# Make directory on EOS and copy tarball from LPC
xrdfs root://cmseos.fnal.gov/ mkdir /store/user/smurthy/condor
xrdcp -f ${localPathToTarball}/${CMSSW}.tgz root://cmseos.fnal.gov//store/user/smurthy/condor/${CMSSW}.tgz
cd -