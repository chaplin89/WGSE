@echo off
REM
REM WGS Extract Reference Library start script (Microsoft Windows)
REM Copyright (C) 2022-2024 Randolph Harr
REM
REM License: GNU General Public License v3 or later
REM A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

REM  Need to be in the installation as CWD could be anywhere if called from the GUI or Symlink
TITLE WGSE Library

set "wgse_FP=%~dp0%"
cd /d "%wgse_FP%"

set "cygbin=%wgse_FP%\cygwin64\bin"

REM Set local PATH path variable for immediate use; registry and PATH not re-read; cannot use variable in PATH
if ";%PATH:cygbin=%;" == ";%PATH%;" ( set "PATH=%cygbin%;%PATH%" )

"%cygbin%\bash.exe" "scripts\library_comman.sh" "%wgse_FP%"
REM "%CYGBIN%\bash.exe" "%~dp0%scripts\get_refgenomes.sh" debug library
pause
