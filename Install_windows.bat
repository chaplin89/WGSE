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

REM Cygwin has Python and Java. Maybe an easier way to install, maintain and use (Unix PATH variable).

REM  Need to be in the installation as CWD could be anywhere
set "WGSEFIN=%~dp0%"
cd /d "%WGSEFIN%"

echo ======================================================================================================
echo WGS Extract v4 Installer for Microsoft Windows Systems
echo(
echo   Will install everything in this installation directory (Python, Java, WGS Extract program, etc)
echo   Bioinformatic tools are natively compiled and available in cygwin64/usr/local/bin
echo(

REM Default URL for the latest release in an online JSON file if local release.json file does not exist
set "track=Alpha"
set "common_url=https://raw.githubusercontent.com/WGSExtract/WGSExtract-Dev/master/"
set "cygwin64version_url=%common_url%/latest-release-%track%.json"
REM Other four URL variables are set in zcommon.sh

setlocal enabledelayedexpansion
set "curlx=%windir%\SysWOW64\curl.exe -fkL"

set "zipdir=cygwin64"                 &  REM Directory stored in the toplevel of the ZIP file
set "reldir=cygwin64"                 &  REM Directory wanting to put the Cygwin64 release into on the local machine
set "latestZIP=cygwin64.zip"          &  REM Name to give ZIP file after download via cURL
set "CYGBIN=%WGSEFIN%cygwin64\bin"    &  REM Full path to CygWin64 bin; once installed
echo CYGBIN %CYGBIN%

REM Cygwin64 Installer and Repository Mirror (albeit repository stored in .zip file)
set "cygsetup=setup-x86_64.exe"
set "cygsetupURL=https://www.cygwin.com/%cygsetup%"
set "mirror=https://cygwin.mirror.constant.com/"    & REM Removed reliance on online mirror; use dir mirror now
set "mirenc=https%%3a%%2f%%2fcygwin.mirror.constant.com%%2f"  & REM Rename to "mirror" after extracting from ZIP

REM todo check http code https://stackoverflow.com/questions/59741113/cmd-batch-if-condition-for-curl-http-response-codes

REM Will add cygwin64\bin to supplied directory then prepend to PATH; sets variable WGSEBIN -- runs elevated
:: call set_WGSEpath.bat "%~dp0"

REM Set local PATH path variable for immediate use; registry and PATH not re-read; cannot use variable in PATH
set "PATH=%CYGBIN%;%PATH%"
REM Possible PATH variable additions needed: python, python\scripts, jre\bin, /usr/local/bin
REM moved before setlocal/endlocal pair; was previously set after installation and befre bash call

REM See https://stackoverflow.com/questions/16042339/os-name-variable
:: for /f "usebackq tokens=1,2 delims==|" %%I in (`wmic os get os_name /format:list`) do 2>NUL set "%%I=%%J"
REM was previously parameter passed to bash; now just use dummy placeholder and calculate inside stage2 script

echo ======================================================================================================
echo Getting latest WGS Extract Cygwin64 package release info
echo(

REM Getting jq.exe simply now before we get it more permanently with Cygwin later
set "jqx=jq.exe"
set "jqURL=https://github.com/stedolan/jq/releases/download/jq-1.6/jq-win64.exe"
if not exist %jqx% (
  echo Downloading jq for Windows ^(json processor^; bootstrap^)
  %curlx% -o %jqx% %jqURL%
)

REM override default online location with local setting; should always exist as delivered with this installer
setlocal enabledelayedexpansion
if exist release.json (
  %jqx% -r .^"release^".^"track^" release.json > temp.txt & set /p track=<temp.txt  & del temp.txt
  %jqx% -r .^"release^".^"!track!URL^" release.json > temp.txt & set /p cygwin64version_url=<temp.txt  & del temp.txt
)

REM Get latest release info to compare against and use URL to download new release; if needed
echo Downloading the latest-release JSON file for track !track!
%curlx% -o latest.json !cygwin64version_url!
if exist latest.json (
  echo Retrieving the latest release version info for package cygwin64 in WGS Extract
  %jqx% -r .^"cygwin64^".^"version^" latest.json > temp.txt  &  set /p  latestVer=<temp.txt & del temp.txt
  %jqx% -r .^"cygwin64^".^"date^"    latest.json > temp.txt  &  set /p latestDate=<temp.txt & del temp.txt
  %jqx% -r .^"cygwin64^".^"URL^"     latest.json > temp.txt  &  set /p  latestURL=<temp.txt & del temp.txt
  echo Found cygwin64 latest release version !latestVer!, date !latestDate!
) else (
  echo Error downloading the latest release JSON file
)

REM del latest.json   & REM leave latest.json for stage2 and common install scripts

REM Check if we need to load the/new CYGWIN64 base package in Windows to bootstrap the rest of the install in BASH
set currentVer=1

if exist cygwin64\ (
  setlocal enabledelayedexpansion

  if exist cygwin64\cygwin64.json (
    %jqx% -r .^"cygwin64^".^"version^" cygwin64\cygwin64.json > temp.txt & set /p  currentVer=<temp.txt & del temp.txt
    %jqx% -r .^"cygwin64^".^"date^"    cygwin64\cygwin64.json > temp.txt & set /p currentDate=<temp.txt & del temp.txt
    %jqx% -r .^"cygwin64^".^"URL^"     cygwin64\cygwin64.json > temp.txt & set /p  currentURL=<temp.txt & del temp.txt
    echo Found cygwin64 current installed version !currentVer!, date !currentDate!

    if !currentVer! GEQ %latestVer% (
      echo Cygwin64 Base v%latestVer% already exists and is up to date ...
      endlocal
      goto postinstall
    )
  )
  endlocal
  set "begin_mesg=Replacing the Cygwin64 Base environment ..."
  set "end_mesg=WGSE Cygwin64 v%latestVer% base installed"
  rmdir cygwin64\ /s/q
) else (
  set "begin_mesg=Installing a new CygWin64 Base v%latestVer% environment (30 second delay for virus check) ..."
  set "end_mesg=WGSE Cygwin64 v%latestVer% base newly installed"
)
echo(

REM Cygwin64 Install / Upgrade routine (just fall into; not called)
:install
echo ======================================================================================================
echo === %begin_mesg%

echo Downloading the Cygwin64 installer %cygsetup%
%curlx% -o %cygsetup% %cygsetupURL%

REM Exclamation in MS Onedrive URL causes issue if delayedexpansion in on; cannot escape it
setlocal disabledelayedexpansion
echo Downloading the Cygwin64 v%latestVer% Install Package (WGS Extract specific; 45 second virus check)
%curlx% -o %latestZIP% "%latestURL%"
endlocal

REM echo cygsetup %cygsetup% %latestZIP% %latestURL%
if exist %cygsetup% (
  if exist %latestZIP% (
    powershell Expand-Archive -LiteralPath "%latestZIP%" -DestinationPath "." -force

    REM The ZIP hard codes %zipdir% inside itself.  Need to rename if %reldir% is different
    if not "%zipdir%" == "%reldir%" (rename %zipdir% %reldir%)
    REM We want to take the http coded mirror directory and simply change it to "mirror"; so no online source checked
    if exist %mirenc%\ (rename %mirenc% mirror)

    echo Starting the Cygwin64 Base Setup -- takes 10 minutes -- see cygwin64\setup.log for a detailed log
    %cygsetup% --root %reldir% --site mirror --only-site --quiet-mode --no-shortcuts ^
      --no-admin --local-package-dir %WGSEFIN% --local-install --categories base ^
      --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,^
libncurses-devel,libcurl-devel,libssl-devel > %reldir%\setup.log

    %CYGBIN%\cygstart.exe /bin/ln.exe -s /cygdrive /mnt
    type %reldir%\etc.skel.bashrc >> %reldir%\etc\skel\.bashrc

    rmdir mirror\ /s/q              & REM Remove the pre-install packages
    rmdir %reldir%\usr\local /s/q   & REM Empty directories there; should be missing for bioinfo version check
    del %cygsetup% %latestZIP%

    echo === %end_mesg%
  ) else (
    echo "*** Failed to download the %latestZIP% archive"
  )
) else (
  echo "*** Failed to download the %cygsetup% executable"
)

:postinstall
REM jq.exe is installed with the cygwin64 base package so delete our temporary local copy
del %jqx%
echo(

if exist "%CYGBIN%\bash.exe" (
  REM Now that we have a working BASH, we install the bioinformatic tools, python, etc in the Stage2 BASH script.
  REM We need to make sure we use our newly availableBASH and not the Win10 supplied one
  echo ======================================================================================================
  echo Starting Stage 2 Windows install in new Cygwin64 BASH environment
  echo(
  "%CYGBIN%\bash.exe" scripts\zinstall_stage2windows.sh dummy

  if %errorlevel% EQU 0 (
    echo ======================================================================================================
    echo(
    echo Congratulations!  You finished installing WGS Extract v4 on MS Windows!
    echo  You can start WGS Extract v4 by clicking the WGSExtract.bat file. Make a softlink, rename it to
    echo  WGSExtract, and place it on your desktop for ease in starting the program.
    echo(
  ) else ( if %errorlevel% EQU 10 (
    REM Restarted the install script due to an upgrade; so exit silently from this one
    goto :EOF
  ) else (
    echo ======================================================================================================
    echo(
    echo Appears there was one or more ERRORS during the WGS Extract v4 install on MS Windows.
    echo  Please scroll back through the command log and look for any errors.
    echo(
  ))
) else (
  echo *** Fatal ERROR cygwin64 failed to install
)

echo Scroll up to first review the log before closing the window.
pause
