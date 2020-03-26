#!/bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lofarbeam/lib
export PATH=$PATH:/opt/dp3/bin
cd ${TESTDIR} && DPPP msin=tNDPPP-generic.MS msout=. steps=[predict] predict.usebeammodel=True predict.elementmodel=oskardipole predict.sourcedb=tNDPPP-generic.MS/sky
