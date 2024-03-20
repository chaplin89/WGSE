@echo off
REM WGS Extract Library Start Script for Ubuntu Linux
REM Copyright (C) 2020-2022 Randolph Harr
REM
REM License: GNU General Public License v3 or later
REM A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

TITLE WGS Extract

REM  Need to be in the installation as CWD could be anywhere
set "wgse_FP=%~dp0%"
cd /d "%wgse_FP%"

echo Starting WGS Extract ...

python\python program\wgsextract.py %*
pause
