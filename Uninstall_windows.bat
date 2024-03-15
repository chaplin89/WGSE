@echo off
REM WGS Extract v4 release
REM Copyright (c) 2021-22 Randy Harr
REM
REM Windows 10/11 native installer (initial, stage 1)
REM
REM Install CygWin64 environment (instead of previous bootstrap only).  Cannot install from CygWin64 installer
REM  as version will not match pre-compiled bioinformatic binaries. So we comment out that approach until
REM  we can resolve that issue. Leave most installs till we start BASH environment (stage 2 and beyond). Easier
REM  there and some in-common with other platforms.
REM

echo ======================================================================================================
echo WGS Extract v4 Uninstaller for Microsoft Windows Systems
echo(
echo   Windows uninstall is simple as all dependent programs are in the WGS Extract installation folder.
echo   So we simply have to delete its install directory and maybe the relocated Reference Library.
echo(

REM  Need to be in the installation as CWD could be anywhere
cd /d %~dp0%

set WGSEFIN=%~dp0%
set reldir=cygwin64
set CYGBIN=%WGSEFIN%\%reldir%\bin

echo ======================================================================================================
echo Finding Reference Library (in case moved in application settings)
echo(

REM Make sure JQ is available
set "jqx=%CYGBIN%\jq.exe"
if not exist %jqx% (
  set jqx=jq.exe
  %curlx% -o %jqx% https://github.com/stedolan/jq/releases/download/jq-1.6/jq-win64.exe
)

set newreflib=
if exist "%USERPROFILE%\.wgsextract" (
  %jqx% -r '."reflib.FP"' "%USERPROFILE%\.wgsextract" > temp.txt & set /p newreflib=<temp.txt & del temp.txt
)
set "reflib=%~dp0%reference"
if exist "%newreflib%\\" (
  echo Original Reference Library in %reflib% but using new Reference Library of %newreflib%
  set "reflib=%newreflib%
)

REM Set local PATH path variable for immediate use; registry and PATH not re-read; cannot use variable in PATH
set "PATH=%CYGBIN%;%PATH%"

REM We need to make sure we use our BASH and not the Win10 supplied one
%CYGBIN%\bash.exe scripts\zuninstall_common.sh "%reflib%
if not exist %WGSEFIN%program/program.json (
  if exist %CYGBIN% (
    REM Wait a short pause for delete to complete; then try to delete the remaining bash.exe left hanging
    sleep 2
    rmdir %WGSEFIN% /s/q
  )
)

