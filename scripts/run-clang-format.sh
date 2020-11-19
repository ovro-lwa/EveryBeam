#!/bin/bash

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# run-clang-format.sh: Formats source code in this repo in accordance with .clang-format file.
#
# To hook this script to pre-commit include the line
# "./scripts/run-clang-format.sh" to .git/hooks/pre-commit
# and make sure pre-commit is an executable shell script.

# Disable globbing
set -e -f

#Script configuration for this repo. Adjust it when copying to a different repo.

#The directory that contains the source files, which clang-format should format.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external)

#The extensions of the source files, which clang-format should format.
SOURCE_EXT=(*.cc *.h)

#End script configuration.

# Detect run environment.
if [ -n "$CI" ]; then
  DRYRUN=" (dry run on CI)"
elif [ -n "$GIT_AUTHOR_DATE" ]; then
  DRYRUN=" (dry run in git hook)"
fi

# print in bold-face
echo -e "\e[1mRunning clang-format$DRYRUN...\e[0m"

# Convert SOURCE_EXT into "-name ext1 -o -name ext2 -o name ext3 ..."
FIND_NAMES="-name ${SOURCE_EXT[0]}"
for i in `seq 1 $((${#SOURCE_EXT[*]} - 1))`; do
  FIND_NAMES+=" -o -name ${SOURCE_EXT[$i]}"
done

# Convert EXCLUDE_DIRS into "-path ./dir1 -prune -o -path ./dir2 -prune -o ..."
FIND_EXCLUDES=
for e in ${EXCLUDE_DIRS[*]}; do
  FIND_EXCLUDES+="-path ./$e -prune -o "
done

cd $SOURCE_DIR
FILES=$(find . $FIND_EXCLUDES -type f \( $FIND_NAMES \) -print)

if [ -n "$DRYRUN" ]; then
  # If the xml has no replacement entries, all files are formatted.
  if clang-format -style=file --output-replacements-xml $FILES |
     grep -q "<replacement "; then
    # Print in bold-face red
    echo -e "\e[1m\e[31mAt least one file is not properly formatted!\e[0m"
    echo -e "\e[1m\e[31mRun $0 for formatting all files!\e[0m"
    exit 1;
  else
    # print in bold-face green
    echo -e "\e[1m\e[32mGreat job, all files are properly formatted!\e[0m"
    exit 0;
  fi
else
  clang-format -i -style=file $FILES
  # print in bold-face
  echo -e "\e[1mSuccessfully formatted all files.\e[0m"
fi