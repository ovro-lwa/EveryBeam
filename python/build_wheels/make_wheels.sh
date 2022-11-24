#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-only

# Script to make python wheels for several versions

SCRIPT_DIR=$(cd $(dirname $0) && pwd)
ROOT_DIR=$(git rev-parse --show-toplevel)

set -euo pipefail
for py_version in 310 39 38 37 36; do
    docker build \
        --build-arg PY_VERSION=${py_version} \
        --file ${SCRIPT_DIR}/py_wheel.docker \
        --tag everybeam-py${py_version} \
        ${ROOT_DIR}
    dockerid=$(docker create everybeam-py${py_version})
    docker cp ${dockerid}:/dist ${ROOT_DIR}
    docker rm ${dockerid}
done
