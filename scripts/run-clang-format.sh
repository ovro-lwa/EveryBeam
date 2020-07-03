#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH
# Move up to parent folder which contains the source
cd ..

# Find all cpp headers ("*.h") and code files ("*.cc") source tree and execute clang-format. Exclude ./external director
find . -path ./external -prune -o -type f \( -name "*.h" -o -name "*.cc" \) -exec clang-format -i -style=file \{\} +
