@echo off
REM
REM WGS Extract Library Start Script for Microsoft Windows
REM Copyright (C) 2022 Randolph Harr
REM
REM License: GNU General Public License v3 or later
REM A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

REM  Need to be in the installation as CWD could be anywhere if called from the GUI or Symlink
cd /d %~dp0%
set "CYGBIN=%~dp0%\cygwin64\bin"

"%CYGBIN%\bash.exe" "%~dp0%\scripts\zlibrary_common.sh" dummy
pause