@echo off
REM WGS Extract Library Start Script for Ubuntu Linux
REM Copyright (C) 2020-2022 Randolph Harr
REM
REM License: GNU General Public License v3 or later
REM A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

REM  Need to be in the installation as CWD could be anywhere if called from the GUI or Symlink
cd /d %~dp0%

echo Starting WGS Extract

python\python program\wgsextract.py
pause