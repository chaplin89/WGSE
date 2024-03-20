#!/usr/bin/env bash
# WGS Extract Program Start Script for Apple MacOS
# Copyright (C) 2020-2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# Find the installation directory (normally part of zcommon.sh but when run standalone, cannot find that script)
WGSEDIR=$(/usr/bin/dirname "${BASH_SOURCE[0]}")   # Get the script location to determine install directory
WGSEABS=$(cd "$WGSEDIR"; pwd -P)                  # By cd'ing to it, resolve any aliases and symlinks
[[ $(/usr/bin/basename "${WGSEABS}") == "scripts" ]] && WGSEFIN=$(/usr/bin/dirname "${WGSEABS}") || WGSEFIN="${WGSEABS}"

cd "${WGSEFIN}"

echo 'Starting WGSExtract...'

/usr/local/bin/python3 "${WGSEFIN}/program/wgsextract.py" "$@"
read -n1 -r -p 'Press enter / return to close this window (first scroll up to review if desired) ...'
