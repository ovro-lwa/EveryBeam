#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Following variables are assumed to be specified in CMakeLists file
# - $LIB_DIR path to pyeverybeam shared library
# - $DATA_DIR path to the (mock) test data
# - $SCRIPTS_DIR path to the (download) scripts

export PYTHONPATH=$LIB_DIR

echo $DATA_DIR
if [ -d $DATA_DIR/LOFAR_LBA_MOCK.ms ]
then
    echo "LBA mock ms exists"
else
    echo "Download LBA mock ms"
    $SCRIPTS_DIR/download_lofar_lba.sh
fi

pytest -s --exitfirst ${DIR}/test_load.py 