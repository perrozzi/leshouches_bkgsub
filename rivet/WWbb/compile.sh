#!/bin/bash
cp WWbb.cc /afs/cern.ch/work/p/perrozzi/private/git/CMSSW_7_5_0_pre5/src/GeneratorInterface/RivetInterface/plugins/
cd /afs/cern.ch/work/p/perrozzi/private/git/CMSSW_7_5_0_pre5/src
eval `scramv1 runtime -sh`
scram b -j16
# export RIVET_ANALYSIS_PATH=${PWD}
cd -

# source /afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67b/rivet/2.2.1/x86_64-slc6-gcc47-opt/rivetenv-genser.sh

# /afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67b/rivet/2.2.1/x86_64-slc6-gcc47-opt/rivetenv.sh

# rivet-buildplugin RivetWWbb.so WWbb.cc

rivet-buildplugin RivetMC_WWBB.so MC_WWBB.cc
rivet --analysis=MC_WWBB ../../cards/HepMC.tt.LO.hepmc2g
