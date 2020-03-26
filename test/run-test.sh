#!/bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lofarbeam/lib
export PATH=$PATH:/opt/dp3/bin
cd ${TESTDIR} && tar xf ${INSTALLDIR}/dp3/DP3-${DP3_VERSION}/DPPP/test/tNDPPP-generic.in_MS.tgz
cd ${TESTDIR} && DPPP msin=tNDPPP-generic.MS msout=. steps=[predict] predict.usebeammodel=True predict.elementmodel=oskardipole predict.sourcedb=tNDPPP-generic.MS/sky
