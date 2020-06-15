#!/bin/sh
# Do check on clang-format, script is inspired on:
# - https://dev.to/10xlearner/formatting-cpp-c-javascript-and-other-stuff-2pof
# - https://gitlab.freedesktop.org/monado/monado/-/blob/master/scripts/format-and-spellcheck.sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl
#
# To hook this script to pre-commit include the line "./scripts/clang-format-check.sh" to .git/hooks/pre-commit
# and make sure pre-commit is an executable shell script.

set -e 
PATCH_NAME=clang-fixes.diff

echo "Running clang-format..."
echo

# Run clang-format from run-clang-format.sh
SCRIPT_PATH=$(dirname "$0")
$SCRIPT_PATH/run-clang-format.sh

# Can't use tee because it hides the exit code
if git diff --patch --exit-code > $PATCH_NAME; then
    echo
    # print in bold-face green
    echo "\e[1m\e[32mclang-format and codespell changed nothing."
    rm -f $PATCH_NAME
    exit 0;
else
    echo
    # Print in bold-face red
    echo "\e[1m\e[31mclang-format made at least one change. Run clang-format before pushing!"
    rm -f $PATCH_NAME
    exit 1;
fi

