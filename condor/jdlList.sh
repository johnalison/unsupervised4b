source makeTarball.sh
mkdir stderr
mkdir stdout
mkdir log
condor_submit wkdtDataSelf.jdl
# condor_submit wkdtDataTT.jdl