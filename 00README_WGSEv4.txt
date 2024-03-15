WGS Extract Beta v4
https://wgsextract.github.io/

Congratulations on downloading the WGS Extract v4 Installer and Program.

There is much more detailed information in the user manual available at
the above link.

Click on the appropriate installer script for your computer :
    Microsoft Win10/11      	Install_windows.bat
    Apple MacOS X/11/12     	Install_macos.command
    Ubuntu Linux (18/120/22)	Install_ubuntu.sh
On Windows or MacOS, Right-Click or CTRL-Click the first time to go
 though the approval process used on your OS. After that, you can simply
 click the script file directly to start. Some OS's have this extra level
 of protection applied to applications not sourced their app store.

If on a MacOS, drag the WGSExtractv4/ folder (result after unzipping) into
 /Applications or maybe your Home folder.  On Windows, drag to
 "C:\Program Files" or your Users Home folder. In both cases, the Home is one
 up from the Downloads folder. You cannot leave the WGS Extract installation
 in the Downloads/ or Documents/ folders; programs cannot work from there.

If upgrading, copy the files and folders from the expanded WGSExtractv4/ folder
 in the same directory as your current installation.  Overwriting or replacing
 any same-named files there.  Then simply (re)run the Installer named above.
 You can rename the installation directory, if desired.

Installers are re-entrant. They will update tools, libraries and packages of
 the WGS Extract program if run again.

We try to make sure a script starts in a visible Terminal / cmd.exe window when
 you click to start it in your GUI / Finder / Explorer.  But it does not always
 work. On all platforms, starting the script from a Terminal / cmd.exe shell
 can assure you get the Terminal window and see any log of actions.  To start
 a script in Win10/11, you "call" a .bat file.  In all others, simply type
 "bash" and then give the full script file name.

If doing a fresh install or major upgrade, the Library* command is called at
 the end of the installation.  You can simply answer "(1) Exit". And then just
 plan to run the Library* command at a later time. If using a CRAM file to start,
 you will need its reference genome immediately. Download the one needed using
 the Library* command; if you know which. Otherwise, the program will mention
 what it cannot find when it is needed. And you can do it then.

Program installation is only 10-15 minutes whereas the Reference Genome Library
 may take 30 minutes to 6 hours to download and configure; depending on your
 internet connection speed, your computer horsepower and how many reference
 genomes you decide to install.

Once installed successfully, start the WGS Extract program by clicking on the
 Start script for your computer (only one of them should exist):
    Microsoft Win10/11      WGSExtract.bat
    Apple MacOS X/11/12     WGSExtract.command
    Ubuntu Linux (LTS)      WGSExtract.sh

Make an alias / shortcut of the Start script file and put it wherever convenient
 (desktop, Application folder, Start menu, etc). How you do this is different
 for each OS. Search the internet or see the Installation chapter in the WGS
 Extract v4 Users Manual for more details: https://bit.ly/3JCyZNa

There is a corresponding Uninstall_xxxx.xxx script for each OS type.

Full source code and documentation at https://wgsextract.github.io/