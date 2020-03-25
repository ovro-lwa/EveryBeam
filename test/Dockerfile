#
# base
#
FROM ubuntu:18.04
RUN apt-get update
RUN apt-get upgrade -y

#
# install astronomy packages
#
RUN apt-get -y install casacore-dev
RUN apt-get -y install python-casacore
RUN apt-get -y install aoflagger-dev
RUN apt-get -y install libcfitsio-dev
RUN apt-get -y install wcslib-dev

#
# install misc packages
#
RUN apt-get -y install wget git make cmake g++

#
# install dependencies
#
RUN apt-get -y install libboost-all-dev
RUN apt-get -y install libhdf5-dev
RUN apt-get -y install libfftw3-dev
RUN apt-get -y install libblas-dev
RUN apt-get -y install liblapack-dev
RUN apt-get -y install libgsl-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install libpng-dev
RUN apt-get -y install libgtkmm-3.0-dev

#
# install packages from source
#
ENV INSTALLDIR /opt

#
# install-lofarbeam
#
RUN mkdir -p ${INSTALLDIR}/lofarbeam/build
RUN cd ${INSTALLDIR}/lofarbeam && git clone https://github.com/lofar-astron/LOFARBeam.git lofarbeam
RUN cd ${INSTALLDIR}/lofarbeam/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/lofarbeam ../lofarbeam
RUN cd ${INSTALLDIR}/lofarbeam/build && make -j
RUN cd ${INSTALLDIR}/lofarbeam/build && make install

#
# install-wsclean
#
RUN mkdir -p ${INSTALLDIR}/wsclean/build
RUN cd ${INSTALLDIR}/wsclean && git clone https://git.code.sf.net/p/wsclean/code wsclean
RUN cd ${INSTALLDIR}/wsclean/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/wsclean -DLOFAR_STATION_RESPONSE_LIB=${INSTALLDIR}/lofarbeam/lib/libstationresponse.so -DLOFAR_STATION_RESPONSE_INCLUDE_DIR=${INSTALLDIR}/lofarbeam/include ${INSTALLDIR}/wsclean/wsclean/wsclean
RUN cd ${INSTALLDIR}/wsclean/build && make -j
RUN cd ${INSTALLDIR}/wsclean/build && make install

#
# install-dp3
#
ENV DP3_VERSION 4.1
RUN mkdir -p ${INSTALLDIR}/dp3/build
RUN cd ${INSTALLDIR}/dp3 && wget --quiet https://github.com/lofar-astron/DP3/archive/v${DP3_VERSION}.tar.gz -O dp3-${DP3_VERSION}.tar.gz
RUN cd ${INSTALLDIR}/dp3 && tar xvf dp3-${DP3_VERSION}.tar.gz
RUN cd ${INSTALLDIR}/dp3/build && cmake ../DP3-${DP3_VERSION} -DCASACORE_ROOT_DIR=${INSTALLDIR}/casacore -DCMAKE_PREFIX_PATH=${INSTALLDIR}/lofarbeam -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/dp3
RUN cd ${INSTALLDIR}/dp3/build && make -j
RUN cd ${INSTALLDIR}/dp3/build && make install

#
# run-test
#
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/lofarbeam/lib
ENV PATH $PATH:/opt/dp3/bin
ENV TESTDIR /test
RUN mkdir -p ${TESTDIR}
RUN cd ${TESTDIR} && tar xf /opt/dp3/DP3-4.1/DPPP/test/tNDPPP-generic.in_MS.tgz
RUN cd ${TESTDIR} && DPPP msin=tNDPPP-generic.MS msout=. steps=[predict] predict.usebeammodel=True predict.elementmodel=oskardipole predict.sourcedb=tNDPPP-generic.MS/sky
