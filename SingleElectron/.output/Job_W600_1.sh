#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
#cd /cms/ekoenig4/MonoZprimeJet/CMSSW_10_2_10/src
#cmsenv
cd ${_CONDOR_SCRATCH_DIR}
./analyze ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8}
