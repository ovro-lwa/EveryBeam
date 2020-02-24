RUN FROM kernsuite/base:5

RUN apt update
RUN apt install -y git build-essential vim casacore-dev cmake python3-distutils libblas-dev liblapack-dev libhdf5-dev
RUN cd
RUN mkdir opt && cd opt
RUN git clone https://git.astron.nl/RD/LOFARBeam
RUN cd LOFARBeam/CMake && git checkout oskar && cd ..
RUN mkdir build && cd build
RUN cmake ..
RUN make -j install

RUN apt install -y libcfitsio-dev libpython-dev libxml2-dev libboost-all-dev libpng-dev aoflagger-dev libgtkmm-3.0-dev
RUN cd ~/opt
RUN git clone https://github.com/lofar-astron/DP3
RUN cd DP3 && git checkout oskar
RUN mkdir build && cd build
RUN cmake -DPYTHON_EXECUTABLE=/usr/bin/python2 ..
RUN make -j install

RUN cd && mkdir oskarbeamtest && cd oskarbeamtest
RUN tar xf ~/opt/DP3/DPPP/test/tNDPPP-generic.in_MS.tgz
RUN DPPP msin=tNDPPP-generic.MS msout=. steps=[predict] predict.usebeammodel=True predict.elementmodel=oskardipole predict.sourcedb=tNDPPP-generic.MS/sky
