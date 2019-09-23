# L1CaloTrigger Phase 2 Jet Algorithm for CMSSW

## Step-1: CMSSW_10_5_0_pre1
````
cmsrel CMSSW_10_5_0_pre1
cd CMSSW_10_5_0_pre1/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_5_0_pre1
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.17.15.1
git cms-addpkg L1Trigger/L1TCommon

scram b -j 8
````
More information can be found here
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_10_5_0_pre1

## Step-2: Tyler's Jet Algorithm
````
git remote add truggles git@github.com:truggles/cmssw.git
git fetch truggles phase2-l1taus_10_5_Xv2
git cms-merge-topic -u truggles:phase2-l1taus_10_5_X_Apr26
````
More information can be found here
https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2L1CaloJetsAndTaus#Central_Integration_Branch

## Step-3: Phase 2 Jet Algorithm Code
````
cd $CMSSW_BASE/src/L1Trigger/
rm -rf L1CaloTrigger/
git clone https://github.com/ekoenig4/L1CaloTrigger.git
````

## Step-4: Crab Jobs
Follow instructions for Tylers code
https://github.com/truggles/L1EGRateStudies/tree/10_5_X_Taus
