# Pi0Tuplizer
Simple photon tuplizer with GEN input
=============================
-----------------------------
setup
-----------------------------
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
if [ -r CMSSW_7_1_25_patch5/src ] ; then
 echo release CMSSW_7_1_25_patch5 already exists
else
scram p CMSSW CMSSW_7_1_25_patch5
fi
cd CMSSW_7_1_25_patch5/src
eval `scram runtime -sh`
git clone https://github.com/zhangzc11/GenPhotonTuplizer.git GenPhoton/GenPhotonTuplizer
scram b
```
-----------------------------
cmsRun
-----------------------------
```
cd python
cmsRun GenPhotonTuplizer_71X_GEN.py
```

-----------------------------
submit crab jobs
-----------------------------
```
crab submit -c crab.py
```


