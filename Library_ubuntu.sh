#!/usr/bin/env bash
# WGS Extract Library Start Script for Ubuntu Linux
# Copyright (C) 2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# Find the installation directory (normally part of zcommon.sh but when run standalone, cannot find that script)
WGSEDIR=$(/usr/bin/dirname "${BASH_SOURCE[0]}")   # Get the script location to determine install directory
WGSEABS=$(cd "$WGSEDIR"; pwd -P)                  # By cd'ing to it, resolve any aliases and symlinks
[[ $(/usr/bin/basename "${WGSEABS}") == "scripts" ]] && WGSEFIN=$(/usr/bin/dirname "${WGSEABS}") || WGSEFIN="${WGSEABS}"

cd "${WGSEFIN}"

if ! source "${WGSEFIN}/scripts/zxterm_ubuntu.sh" ; then    # Start an xterm if not in a Terminal to start
  echo "*** ERROR sourcing scripts/zxterm_ubuntu.sh.sh in ${BASH_SOURCE[0]}. Exiting."
  exit 1
fi

/usr/bin/env bash "${WGSEFIN}/scripts/get_refgenomes.sh" library
