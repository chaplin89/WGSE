# coding: utf8
#
# Microarray File generation subsystem (formally known as the Autosomal File Format generation)
#
# Part of the
# WGS Extract (https://wgse.bio/) system
# Todo modify to allow stand alone operation but may still be dependent on WGS Extract support utilities, settings
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys      # for argv[] if called as standalone program
import os       # for os.path
# import shutil   # for shutil.move
from functools import partial
from tkinter import Toplevel, BooleanVar, W, Checkbutton, Label, LabelFrame, Frame, LEFT, CENTER
try:
    from tkmacosx import Button         # Only way to get colored background buttons on MacOS
except ImportError:
    from tkinter import Button

from utilities import nativeOS, DEBUG, unquote, wgse_message, offset
from commandprocessor import run_bash_script
# from mainwindow import mainwindow_resume      # Localized inside cancel_miroarray_window() due to loop
import settings as wgse


# Todo Microarray Format Window update -- as tab to main window? New standalone program and class?

"""###################################################################################################################
  Microarray File Generator subsystem / window (aka Autosomal File Formats)
    This original microarray files generator was developed in WGSE Extract by City Farmer.  The largest, most
    complete and unique contribution of the WGS Extract tool.  It grew out of the basic concept in Extract23 that
    took a 30x WGS and generated a 23andMe v3 file.  Except based on a generated AllSNPs template file (in build 37 
    reference format) which is then subsetted for each of vendor formats (using a module "aconv" that is developed 
    here). Also does a liftover (based on wrapping pyliftover rewrite of UCSC's liftover in hg38tohg19.py module) 
    from Build 38 to 37 of the AllSNPs file generated if original BAM/CRAM is in Build 38 format.
"""

minAllZipsize = 2500000         # Minimum All SNP ZIP size; pretty stable at size of template
minAllSNPSize = 5000000         # Parameterized minimum All SNP text size that is considered valid (Bytes)
# Most often is 17.5 MB. Had it set at 5 MB but Teemu wants it smaller for aDNA samples he runs through to get some
#  ad-mixture analysis from. Cannot recall why was set at 5 MB; maybe due to some common error we were encountering.
# Ditto for AllSNPsSize when not all the values are available to be called.
# Todo include no calls in vendor files so size is more static and determinable

# IMPORTANT: Entries here must match exactly as used in Button definition and aconv.py; order is not important
targets = [
    # all multi- / combined- file targets (vendor, Ancient, WGSE)
    'AllSNPs', 'CombinedKit_v3', 'CombinedKit_v5', '23andMe_V35', '23andMe_SNPs_API', "1240K+HO",
    # Vendor targets (genealogy)
    '23andMe_V3', '23andMe_V4', '23andMe_V5', 'Ancestry_V1', 'Ancestry_V2',
    'FTDNA_V1_Affy', 'FTDNA_V2', 'FTDNA_V3', 'LivDNA_V1', 'LivDNA_V2', 'MyHeritage_V1', 'MyHeritage_V2',
    # Others (health and wellness)
    "MTHFRGen", "TellMeGen", "Genera", "meuDNA", "SelfDecode",
    # Ancient DNA
    "1240K", "HOv1"
]
# todo need to use the above to create dictionaries of buttons and results instead of enumerated below
#  possibly adding to the above attributes such as recommended, file extension type, template name, etc
#  so the code is completely independent of any particular instance.  Allowing the above to be updated,
#  along with a new template added, to add a new vendor button. Not sure how to handle GUI selection then.

# Any use of these individual variables is not order important. Just kept the same order for consistency
allSNPsButton       = None   ;   allSNPsButtonResult       = None
combinedv3Button    = None   ;   combinedv3ButtonResult    = None
combinedv5Button    = None   ;   combinedv5ButtonResult    = None

g23andmev3Button    = None   ;   g23andmev3ButtonResult    = None
g23andmev4Button    = None   ;   g23andmev4ButtonResult    = None
g23andmev5Button    = None   ;   g23andmev5ButtonResult    = None
g23andmev35Button   = None   ;   g23andmev35ButtonResult   = None
g23andmeAPIButton   = None   ;   g23andmeAPIButtonResult   = None

ancestryv1Button    = None   ;   ancestryv1ButtonResult    = None
ancestryv2Button    = None   ;   ancestryv2ButtonResult    = None

ftdnav1Button       = None   ;   ftdnav1ButtonResult       = None
ftdnav2Button       = None   ;   ftdnav2ButtonResult       = None
ftdnav3Button       = None   ;   ftdnav3ButtonResult       = None

livdnav1Button      = None   ;   livdnav1ButtonResult      = None
livdnav2Button      = None   ;   livdnav2ButtonResult      = None

myheritagev1Button  = None   ;   myheritagev1ButtonResult  = None
myheritagev2Button  = None   ;   myheritagev2ButtonResult  = None

mthfrGenButton      = None   ;   mthfrGenButtonResult      = None
tellmegenButton     = None   ;   tellmegenButtonResult     = None
generaButton        = None   ;   generaButtonResult        = None
meuDNAButton        = None   ;   meuDNAButtonResult        = None
selfdecodeButton    = None   ;   selfdecodeButtonResult    = None

g1240KButton        = None   ;   g1240KButtonResult        = None
hov1Button          = None   ;   hov1ButtonResult          = None
g1240hoButton       = None   ;   g1240hoButtonResult       = None

checkbuttons_results = []
microarrayWindow = None
generateSelectedFilesButton = None


def button_microarray_window():
    """
        Main entry into the microarray file generator from WGSE (or as a standalone with locally recreated settings).
        Generates main selection window and then waits on User direction. Utilizes selected (in settings) BAM and
        Output directory area.  Final output are microarray files selected by user.
        Note: if finds an existing AllSNPs file, will start with that instead of recreating it
    """
    global allSNPsButton, combinedv3Button, combinedv5Button
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav1Button, ftdnav2Button, ftdnav3Button
    global livdnav1Button, livdnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, tellmegenButton, generaButton, meuDNAButton, selfdecodeButton
    global g1240KButton, hov1Button, g1240hoButton

    global allSNPsButtonResult, combinedv3ButtonResult, combinedv5ButtonResult
    global g23andmev3ButtonResult, g23andmev4ButtonResult, g23andmev5ButtonResult, \
           g23andmeAPIButtonResult, g23andmev35ButtonResult
    global ancestryv1ButtonResult, ancestryv2ButtonResult, ftdnav1ButtonResult, ftdnav2ButtonResult, ftdnav3ButtonResult
    global livdnav1ButtonResult, livdnav2ButtonResult, myheritagev1ButtonResult, myheritagev2ButtonResult
    global mthfrGenButtonResult, tellmegenButtonResult, generaButtonResult, meuDNAButtonResult, selfdecodeButtonResult
    global g1240KButtonResult, hov1ButtonResult, g1240hoButtonResult

    global generateSelectedFilesButton, checkbuttons_results, microarrayWindow

    font = wgse.fonts.table

    # Create main selection window for the Autosomal microarray file generator
    # wgse.window.withdraw()
    microarrayWindow = Toplevel(wgse.window)
    microarrayWindow.transient()
    microarrayWindow.protocol("WM_DELETE_WINDOW", cancel_miroarray_window)
    microarrayWindow.title(wgse.lang.i18n['SelectAutosomalFormatsTitle'])
    microarrayWindow.geometry(wgse.microarr_winsize + offset(wgse.window.geometry(), 50))
    microarrayWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    microarrayWindow.columnconfigure(0, weight=1)
    microarrayWindow.rowconfigure(0, weight=1)

    # Create Frames for each grouping and the row of buttons across the bottom
    titleFrame = Frame(microarrayWindow)
    titleFrame.grid(row=1, padx=5, pady=3, sticky="nw")
    titleFrame.columnconfigure(0, weight=1)
    titleFrame.rowconfigure(0, weight=1)

    combinedFrame = LabelFrame(microarrayWindow, text=wgse.lang.i18n['FrameCombinedKit'], font=font['14'])
    combinedFrame.grid(row=2, padx=5, pady=3, sticky='W')
    combinedFrame.columnconfigure(0, weight=1)
    combinedFrame.rowconfigure(1, weight=1)

    vendorFrame = LabelFrame(microarrayWindow, text=wgse.lang.i18n['FrameVendors'], font=font['14'])
    vendorFrame.grid(row=3, padx=5, pady=3, sticky='W')
    vendorFrame.columnconfigure(0, weight=1)
    vendorFrame.rowconfigure(1, weight=1)

    othersFrame = LabelFrame(microarrayWindow, text=wgse.lang.i18n['FrameOthers'], font=font['14'])
    othersFrame.grid(row=4, padx=5, pady=3, sticky='W')
    othersFrame.columnconfigure(0, weight=1)
    othersFrame.rowconfigure(1, weight=1)
    
    aDNAFrame = LabelFrame(microarrayWindow, text=wgse.lang.i18n['FrameaDNA'], font=font['14'])
    aDNAFrame.grid(row=5, padx=5, pady=3, sticky='W')
    aDNAFrame.columnconfigure(0, weight=1)
    aDNAFrame.rowconfigure(1, weight=1)

    buttonsFrame = Frame(microarrayWindow)
    buttonsFrame.grid(row=6, padx=5, pady=3, sticky='W')
    buttonsFrame.columnconfigure(0, weight=1)
    buttonsFrame.rowconfigure(1, weight=1)

    # Generate select box buttons for each available format in each frame
    checkbuttons_results = []

    #  Title frame --------------------------------------------------------------------------------------------------
    Label(titleFrame, text=wgse.lang.i18n['PleaseSelectAutosomalFormats'], font=font['16b'],
          anchor="nw", justify=CENTER).grid(row=0, padx=5, pady=5)

    Label(titleFrame, text=wgse.BAM.disp_FBS, font=font['12']).grid(row=1, padx=5, pady=0)

    Label(titleFrame, text=wgse.lang.i18n['ExplainAutosomalFormats'], font=font['12'],
          anchor="nw", justify=LEFT).grid(row=2, padx=5, pady=5)

    # General idea is for each frame, first column is white space, 2nd column is vendor / row label, then columns
    #  of check mark boxes to select those variations for that vendor / row title

    #  CombinedKit frame --------------------------------------------------------------------------------------------
    crow = 0
    Label(combinedFrame, text="     ", font=font['14'], anchor=W).grid(
        row=crow, column=0, padx=5, pady=2)   # White space before labels (no other labels in column)

    Label(combinedFrame, text=wgse.lang.i18n['WGSExtract'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    allSNPsButtonResult = BooleanVar()  ;  allSNPsButtonResult.set(True)  ;  checkbuttons_results.append(True)
    allSNPsButton = Checkbutton(combinedFrame, variable=allSNPsButtonResult, font=font['14'], anchor=W,
                                text=wgse.lang.i18n['vall'], command=partial(click_checkbox, "AllSNPs"))
    allSNPsButton.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    combinedv3ButtonResult = BooleanVar()  ;  combinedv3ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    combinedv3Button = Checkbutton(combinedFrame, variable=combinedv3ButtonResult, font=font['14'], anchor=W,
                                   text=wgse.lang.i18n['v3'] + '*', command=partial(click_checkbox, "CombinedKit_v3"))
    combinedv3Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W')

    combinedv5ButtonResult = BooleanVar()  ;  combinedv5ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    combinedv5Button = Checkbutton(combinedFrame, variable=combinedv5ButtonResult, font=font['14'], anchor=W,
                                   text=wgse.lang.i18n['v5'], command=partial(click_checkbox, "CombinedKit_v5"))
    combinedv5Button.grid(row=crow, column=4, padx=1, pady=1, sticky='W'); crow += 1

    Label(combinedFrame, text=wgse.lang.i18n['23andMe'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    g23andmev35ButtonResult = BooleanVar()  ;  g23andmev35ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g23andmev35Button = Checkbutton(combinedFrame, variable=g23andmev35ButtonResult, font=font['14'], anchor=W,
                                    text=f"{wgse.lang.i18n['v3']} + {wgse.lang.i18n['v5']} *",
                                    command=partial(click_checkbox, "23andMe_V35"))
    g23andmev35Button.grid(row=crow, column=2, columnspan=2, padx=1, pady=1, sticky='W')

    g23andmeAPIButtonResult = BooleanVar()  ;  g23andmeAPIButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g23andmeAPIButton = Checkbutton(combinedFrame, variable=g23andmeAPIButtonResult, font=font['14'], anchor=W,
                                    text=wgse.lang.i18n['23andMeFutureSNPs'],
                                    command=partial(click_checkbox, "23andMe_SNPs_API"))
    g23andmeAPIButton.grid(row=crow, column=4, columnspan=2, padx=1, pady=1, sticky='W'); crow += 1

    Label(combinedFrame, text=wgse.lang.i18n['ReichLab'] + '           ', font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    g1240hoButtonResult = BooleanVar()  ;  g1240hoButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g1240hoButton = Checkbutton(combinedFrame, variable=g1240hoButtonResult, font=font['14'], anchor=W,
                                text=wgse.lang.i18n['aDNA1240K+HO'],
                                command=partial(click_checkbox, "1240K+HO"))
    g1240hoButton.grid(row=crow, column=2, columnspan=4, padx=1, pady=1, sticky='W'); crow += 1

    # Vendor frame --------------------------------------------------------------------------------------------------
    crow = 0
    Label(vendorFrame, text="     ", font=font['14'], anchor=W).grid(
        row=crow, column=0, padx=5, pady=2)   # White space before labels (no other labels in column)

    # 23andMe versions (2 rows)
    Label(vendorFrame, text=wgse.lang.i18n['23andMe'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    g23andmev3ButtonResult = BooleanVar()  ;  g23andmev3ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g23andmev3Button = Checkbutton(vendorFrame, variable=g23andmev3ButtonResult, font=font['14'], anchor=W,
                                   text=wgse.lang.i18n['v3'] + ' *', command=partial(click_checkbox, "23andMe_V3"))
    g23andmev3Button.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    g23andmev4ButtonResult = BooleanVar()  ;  g23andmev4ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g23andmev4Button = Checkbutton(vendorFrame, variable=g23andmev4ButtonResult, font=font['14'], anchor=W,
                                   text=wgse.lang.i18n['v4'], command=partial(click_checkbox, "23andMe_V4"))
    g23andmev4Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W')

    g23andmev5ButtonResult = BooleanVar()  ;  g23andmev5ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g23andmev5Button = Checkbutton(vendorFrame, variable=g23andmev5ButtonResult, font=font['14'], anchor=W,
                                   text=wgse.lang.i18n['v5'] + ' *', command=partial(click_checkbox, "23andMe_V5"))
    g23andmev5Button.grid(row=crow, column=4, padx=1, pady=1, sticky='W'); crow += 1

    # AncestryDNA
    Label(vendorFrame, text=wgse.lang.i18n['AncestryDNA'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    ancestryv1ButtonResult = BooleanVar()  ;  ancestryv1ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    ancestryv1Button = Checkbutton(vendorFrame, variable=ancestryv1ButtonResult, font=font['14'], anchor=W,
                                   text='v1',
                                   command=partial(click_checkbox, "Ancestry_V1"))
    ancestryv1Button.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    ancestryv2ButtonResult = BooleanVar()  ;  ancestryv2ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    ancestryv2Button = Checkbutton(vendorFrame, variable=ancestryv2ButtonResult, font=font['14'], anchor=W,
                                   text='v2',
                                   command=partial(click_checkbox, "Ancestry_V2"))
    ancestryv2Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W'); crow += 1

    # FamilyTreeDNA
    Label(vendorFrame, text=wgse.lang.i18n['FamilyTreeDNA'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    ftdnav1ButtonResult = BooleanVar()  ;  ftdnav1ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    ftdnav1Button = Checkbutton(vendorFrame, variable=ftdnav1ButtonResult, font=font['14'], anchor=W,
                                text='v1',
                                command=partial(click_checkbox, "FTDNA_V1_Affy"))
    ftdnav1Button.grid(row=crow, column=2, padx=1, pady=1, sticky='W')
    
    ftdnav2ButtonResult = BooleanVar()  ;  ftdnav2ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    ftdnav2Button = Checkbutton(vendorFrame, variable=ftdnav2ButtonResult, font=font['14'], anchor=W,
                                text='v2',
                                command=partial(click_checkbox, "FTDNA_V2"))
    ftdnav2Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W')

    ftdnav3ButtonResult = BooleanVar()  ;  ftdnav3ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    ftdnav3Button = Checkbutton(vendorFrame, variable=ftdnav3ButtonResult, font=font['14'], anchor=W,
                                text='v3',
                                command=partial(click_checkbox, "FTDNA_V3"))
    ftdnav3Button.grid(row=crow, column=4, padx=1, pady=1, sticky='W'); crow += 1

    # LivingDNA
    Label(vendorFrame, text=wgse.lang.i18n['LivingDNA'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    livdnav1ButtonResult = BooleanVar()  ;  livdnav1ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    livdnav1Button = Checkbutton(vendorFrame, variable=livdnav1ButtonResult, font=font['14'], anchor=W,
                                 text='v1',
                                 command=partial(click_checkbox, "LivDNA_V1"))
    livdnav1Button.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    livdnav2ButtonResult = BooleanVar()  ;  livdnav2ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    livdnav2Button = Checkbutton(vendorFrame, variable=livdnav2ButtonResult, font=font['14'], anchor=W,
                                 text='v2',
                                 command=partial(click_checkbox, "LivDNA_V2"))
    livdnav2Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W'); crow += 1

    # MyHeritage
    Label(vendorFrame, text=wgse.lang.i18n['MyHeritage'], font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    myheritagev1ButtonResult = BooleanVar()  ;  myheritagev1ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    myheritagev1Button = Checkbutton(vendorFrame, variable=myheritagev1ButtonResult, font=font['14'], anchor=W,
                                     text='v1',
                                     command=partial(click_checkbox, "MyHeritage_V1"))
    myheritagev1Button.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    myheritagev2ButtonResult = BooleanVar()  ;  myheritagev2ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    myheritagev2Button = Checkbutton(vendorFrame, variable=myheritagev2ButtonResult, font=font['14'], anchor=W,
                                     text='v2',
                                     command=partial(click_checkbox, "MyHeritage_V2"))
    myheritagev2Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W'); crow += 1

    # Others frame ---------------------------------------------------------------------------------------
    crow = 0
    Label(othersFrame, text="                                 ", font=font['14'], anchor=W).grid(
        row=crow, column=0, padx=5, pady=2)   # White space before checkboxes (no other labels in column)

    # TellMeGen
    tellmegenButtonResult = BooleanVar()  ;  tellmegenButtonResult.set(False)  ;  checkbuttons_results.append(False)
    tellmegenButton = Checkbutton(othersFrame, variable=tellmegenButtonResult, font=font['14'], anchor=W,
                                  text='TellMeGen (ES)',
                                  command=partial(click_checkbox, "TellMeGen"))
    tellmegenButton.grid(row=crow, column=1, padx=2, pady=1, sticky='W')

    # SelfDecode
    selfdecodeButtonResult = BooleanVar()  ; selfdecodeButtonResult.set(False)  ;  checkbuttons_results.append(False)
    selfdecodeButton = Checkbutton(othersFrame, variable=selfdecodeButtonResult, font=font['14'], anchor=W,
                                   text='SelfDecode (US)',
                                   command=partial(click_checkbox, "SelfDecode"))
    selfdecodeButton.grid(row=crow, column=2, padx=1, pady=1, sticky='W'); crow += 1

    # Genera
    generaButtonResult = BooleanVar()  ;  generaButtonResult.set(False)  ;  checkbuttons_results.append(False)
    generaButton = Checkbutton(othersFrame, variable=generaButtonResult, font=font['14'], anchor=W,
                               text='Genera (BR)',
                               command=partial(click_checkbox, "Genera"))
    generaButton.grid(row=crow, column=1, padx=1, pady=1, sticky='W')

    # MeuDNA
    meuDNAButtonResult = BooleanVar()  ;  meuDNAButtonResult.set(False)  ;  checkbuttons_results.append(False)
    meuDNAButton = Checkbutton(othersFrame, variable=meuDNAButtonResult, font=font['14'], anchor=W,
                               text='meuDNA (BR)',
                               command=partial(click_checkbox, "meuDNA"))
    meuDNAButton.grid(row=crow, column=2, padx=1, pady=1, sticky='W'); crow += 1

    # MTHFR Genetics
    mthfrGenButtonResult = BooleanVar()  ;  mthfrGenButtonResult.set(False)  ;  checkbuttons_results.append(False)
    mthfrGenButton = Checkbutton(othersFrame, variable=mthfrGenButtonResult, font=font['14'], anchor=W,
                                 text='MTHFR Genetics (UK)',
                                 command=partial(click_checkbox, "MTHFRGen"))
    mthfrGenButton.grid(row=crow, column=1, columnspan=2, padx=1, pady=1, sticky='W')

    # Ancient DNA frame ----------------------------------------------------------------------------------
    crow = 0
    Label(aDNAFrame, text="     ", font=font['14'], anchor="w").grid(
        row=crow, column=0, padx=5, pady=2)   # White space before labels (no other labels in column)

    Label(aDNAFrame, text=wgse.lang.i18n['ReichLab'] + '           ', font=font['14'], anchor=W).grid(
        row=crow, column=1, padx=1, pady=1, sticky='W')

    g1240KButtonResult = BooleanVar()  ;  g1240KButtonResult.set(False)  ;  checkbuttons_results.append(False)
    g1240KButton = Checkbutton(aDNAFrame, variable=g1240KButtonResult, font=font['14'], anchor=W,
                               text=wgse.lang.i18n['aDNA1240K'],
                               command=partial(click_checkbox, "1240K"))
    g1240KButton.grid(row=crow, column=2, padx=1, pady=1, sticky='W')

    hov1ButtonResult = BooleanVar()  ;  hov1ButtonResult.set(False)  ;  checkbuttons_results.append(False)
    hov1Button = Checkbutton(aDNAFrame, variable=hov1ButtonResult, font=font['14'], anchor=W,
                             text=wgse.lang.i18n['aDNAHOv1'],
                             command=partial(click_checkbox, "HOv1"))
    hov1Button.grid(row=crow, column=3, padx=1, pady=1, sticky='W'); crow += 1

    # Generate action buttons across the bottom row
    selectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectEverything'],
                                        font=font['14'], command=button_select_all)
    selectAllFileFormatsButton.grid(column=0, row=0, padx=5, pady=3)

    selectRecFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectRecommended'],
                                        font=font['14'], command=button_select_recommended)
    selectRecFileFormatsButton.grid(column=1, row=0, padx=5, pady=3)

    deSelectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonDeselectEverything'],
                                          font=font['14'], command=button_deselect_all)
    deSelectAllFileFormatsButton.grid(column=2, row=0, padx=5, pady=3)

    generateSelectedFilesButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonGenerateSelectedFiles'],
                                         font=font['14'], command=button_generate_selected, state="disabled")
    generateSelectedFilesButton.grid(column=4, row=0, padx=5, pady=3)

    # Had previously assumed user knew they could click window exit in upper bar; now provide explicit exit button
    exitButton = Button(buttonsFrame, text=wgse.lang.i18n['CloseWindow'], font=font['14'],
                        command=cancel_miroarray_window)
    exitButton.grid(column=5, row=0, padx=5, pady=3)

    # Select Recommended by default and wait for an action button (generate or exit)
    button_select_recommended()
    microarrayWindow.update()
    microarrayWindow.deiconify()    # Update not showing recommended at start after grab_set() added
    microarrayWindow.grab_set()
    microarrayWindow.wait_window()


# noinspection PyUnresolvedReferences
def button_select_recommended():
    """ Routine to setup Recommended selections (default on entry; can be restored by user button).  """
    global checkbuttons_results, generateSelectedFilesButton, targets
    global combinedv3Button, g23andmev3Button, g23andmev5Button, g23andmev35Button

    button_deselect_all()
    combinedv3Button.select()   ;  checkbuttons_results[targets.index("CombinedKit_v3")] = True
    g23andmev3Button.select()   ;  checkbuttons_results[targets.index("23andMe_V3")]  = True
    g23andmev5Button.select()   ;  checkbuttons_results[targets.index("23andMe_V5")]  = True
    g23andmev35Button.select()  ;  checkbuttons_results[targets.index("23andMe_V35")] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_select_all():
    global checkbuttons_results, generateSelectedFilesButton
    global allSNPsButton, combinedv3Button, combinedv5Button
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav1Button, ftdnav2Button, ftdnav3Button
    global livdnav1Button, livdnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, tellmegenButton, generaButton, meuDNAButton, selfdecodeButton
    global g1240KButton, hov1Button, g1240hoButton

    allSNPsButton.select()
    combinedv3Button.select()
    combinedv5Button.select()

    g23andmev3Button.select()
    g23andmev4Button.select()
    g23andmev5Button.select()
    g23andmeAPIButton.select()
    g23andmev35Button.select()
    ancestryv1Button.select()
    ancestryv2Button.select()
    ftdnav1Button.select()
    ftdnav2Button.select()
    ftdnav3Button.select()
    livdnav1Button.select()
    livdnav2Button.select()
    myheritagev1Button.select()
    myheritagev2Button.select()

    mthfrGenButton.select()
    tellmegenButton.select()
    generaButton.select()
    meuDNAButton.select()
    selfdecodeButton.select()

    g1240KButton.select()
    hov1Button.select()
    g1240hoButton.select()

    for i in range(len(checkbuttons_results)):
        checkbuttons_results[i] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_deselect_all():
    global checkbuttons_results, generateSelectedFilesButton
    global allSNPsButton, combinedv3Button, combinedv5Button
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav1Button, ftdnav2Button, ftdnav3Button
    global livdnav1Button, livdnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, tellmegenButton, generaButton, meuDNAButton, selfdecodeButton
    global g1240KButton, hov1Button, g1240hoButton

    # allSNPsButton.deselect()
    combinedv3Button.deselect()
    combinedv5Button.deselect()

    g23andmev3Button.deselect()
    g23andmev4Button.deselect()
    g23andmev5Button.deselect()
    g23andmeAPIButton.deselect()
    g23andmev35Button.deselect()
    ancestryv1Button.deselect()
    ancestryv2Button.deselect()
    ftdnav1Button.deselect()
    ftdnav2Button.deselect()
    ftdnav3Button.deselect()
    livdnav1Button.deselect()
    livdnav2Button.deselect()
    myheritagev1Button.deselect()
    myheritagev2Button.deselect()

    mthfrGenButton.deselect()
    tellmegenButton.deselect()
    generaButton.deselect()
    meuDNAButton.deselect()
    selfdecodeButton.deselect()

    g1240KButton.deselect()
    hov1Button.deselect()
    g1240hoButton.deselect()

    for i in range(len(checkbuttons_results)):
        if i != targets.index("AllSNPs"):
            checkbuttons_results[i] = False
    # generateSelectedFilesButton.configure(state="disabled")   # generate button always enabled as AllSNPs always enabled


def click_checkbox(kit):
    global checkbuttons_results, generateSelectedFilesButton, targets, allSNPsButton

    i = targets.index(kit)
    if kit != "AllSNPs":         # Normal mode; toggle saved state. Button auto selected / deselected by tkinter
        checkbuttons_results[i] = not checkbuttons_results[i]
    else:                        # Leave AllSNPs always checked now so it is never deleted by WGSE
        checkbuttons_results[i] = True
        allSNPsButton.select()

    number_of_selected_buttons = 0
    for i in range(len(checkbuttons_results)):
        if checkbuttons_results[i]:
            number_of_selected_buttons += 1

    state = "normal" if number_of_selected_buttons > 0 else "disabled"
    generateSelectedFilesButton.configure(state=state)


def get_target_type_suffix(target_type_name_all):
    """ Get the file suffix type for the GLOBAL parameter target_type_name_all. """
    # Note: TellMeGen extension is CSV but internally a TSV
    if any(x in target_type_name_all for x in ["FTDNA", "MyHeritage", "TellMeGen", "meuDNA", "Genera"]):
        return ".csv"
    elif any(x in target_type_name_all for x in ["CombinedKit", "23andMe", "Ancestry", "LivDNA", "MTHFR", "SelfDecode", "1240K", "HO"]):
        return ".txt"
    else:
        return ".err"                             # REH 20Mar2020 Refactor


def button_generate_selected():
    """ Run Microarray format generator (originally based on Extract23) to generate selected formats. """
    global microarrayWindow, targets

    if wgse.BAM.chrom_types['A'] == 0:
        # Need at least Autosomes present to create microarray file. Button should not have been available but ...
        # Todo should provide a pop-up warning message that doing nothing here
        cancel_miroarray_window()
        return

    # microarrayWindow.withdraw()

    # Call internal button to reuse or create a new AllSNPs file in the Output Diretory from the current BAM file
    reused = _button_AllSNPs(microarrayWindow)
    if reused is None:
        pass    # Todo add error message and exit

    # By this point, we should have a .txt AllSNPs file (check_valid_microarray has been run and returned True)

    out_FPB = wgse.outdir.FPB
    out_oFPB = wgse.outdir.oFPB
    out_FPB_AllSNPs = f'{wgse.outdir.FPB}_AllSNPs'
    AllSNPsTXT_oFN = nativeOS(f'{out_FPB_AllSNPs}.txt')

    pythonx = wgse.python3x_qFN
    aconv = f'"{wgse.prog_FP}aconv.py"'
    zipx = wgse.zipx_qFN

    # Need AllSNPs VCF template to compare modified times to see if existing vendor file is out of date
    refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, type="Microarray")
    if "error" in refVCFtab_qFN or not os.path.isfile(nativeOS(unquote(refVCFtab_qFN))):
        wgse_message("error", "errRefLibMissingFileTitle", True,
                     wgse.lang.i18n["errRefLibMissingAllTable"].replace("{{FILE}}", refVCFtab_qFN))
        return

    # Now for EACH format selected (other than AllSNPs), add a run of aconv.py to take and subset the
    #  AllSNPs file to create the respective format. Note that we are simply running aconv.py as a stand-alone
    #  python program; not importing and calling it.
    # Todo modify to simply call aconv routine, within same python thread, and use python zip command
    #  speedy enough in python to not need a subprocess / shell call.

    commands = ""
    stop = len(checkbuttons_results) - (0 if wgse.DEBUG_MODE else 7)    # New formats not quite ready for primetime
    for i in range(stop):

        out_oFPB = f'{out_oFPB}_{targets[i]}'        # Only generate and keep compressed form now
        if i != targets.index("AllSNPs") and checkbuttons_results[i] and \
           not _check_valid_microarray(out_oFPB, refVCFtab_qFN):

            suffix = get_target_type_suffix(targets[i])
            out_uncompressed_qFN = f'"{out_FPB}_{targets[i]}{suffix}"'
            out_compressed_qFN = f'"{out_FPB}_{targets[i]}.zip"'
            commands += (
                f'echo "Generating microarray file for format {targets[i]}"\n'
                f'{pythonx} {aconv} {targets[i]} "{out_FPB_AllSNPs}.txt" "{out_FPB}" "{wgse.reflib.cma_FP}"\n'
                f'{zipx} -mj {out_compressed_qFN} {out_uncompressed_qFN}\n'  # Do a move; delete uncompressed
            )

    if commands:    # Could end up having nothing to do; check first
        run_bash_script("ButtonMicroarrayDNA", commands, parent=microarrayWindow)

    else:       # Report whether reused AllSNPs or not
        DEBUG("*** INFO: No microarray vendor format to regenerate; all already exist and are current.")

    # Note: always have checkbuttons_results[AllSNPs] set now so will never delete unless invalid size
    if not (reused or checkbuttons_results[targets.index("AllSNPs")]) or \
        (os.path.isfile(AllSNPsTXT_oFN) and os.path.getsize(AllSNPsTXT_oFN) < minAllSNPSize):
        wgse.tempf.list.append(AllSNPsTXT_oFN)

    # Todo generate a stats window pop-up after running the generator.

    # button_microarray_stats()

    cancel_miroarray_window()


def _check_valid_microarray(oFPB, refVCF_qFN):
    """ Check for a valid microarray file. Needs refVCF to make sure newer than that.  """

    # Decide which file to use to check (either can exist; but choose .txt over .zip
    oFNt = oFPB + '.txt'
    oFNz = oFPB + '.zip'
    if os.path.exists(oFNt) and os.path.isfile(oFNt):
        oFN = oFNt
    elif os.path.exists(oFNz) and os.path.isfile(oFNz):
        oFN = oFNz
    else:
        return False
    oFN_mt = os.path.getmtime(oFN)

    # Determine minimum "valid" size based on whether AllSNPs file or a vendor file
    # todo can be more aggresive (tighter size value) if make microarray files "all call" for each vendor.
    #  aDNA community has very small microarray files that are below any limit we may set here.
    if "_AllSNPs" in oFN:               # AllSNPs file
        checksize = minAllSNPSize

    else:                               # Vendor file
        checksize = minAllZipsize

        # If a vendor file, check now for a modified time as compared to the AllSNPs file
        AllSNPsTXT_oFN = nativeOS(f'{wgse.outdir.FPB}_AllSNPs.txt')  # Only use TXT form for AllSNPs here; never ZIP
        if not os.path.exists(AllSNPsTXT_oFN) or oFN_mt < os.path.getmtime(AllSNPsTXT_oFN):
            return False

    # Return result of final checks (modified time against recVCF template and BAM file; also file size)
    return oFN_mt > os.path.getmtime(nativeOS(unquote(refVCF_qFN))) and \
           oFN_mt > os.path.getmtime(wgse.BAM.file_oFN) and \
           os.path.getsize(oFN) > checksize


def _button_AllSNPs(parent_window, parallel=False):
    """
      AllSNPs merged Microarray File Creator (internal button). If returned None, then error and did not generate and
      likely does not exist. Else returns whether reusing existing file (True) or created a valid new one (False).

      Loosely, originally based on the Extract23 script: https://github.com/tkrahn/extract23/blob/master/extract23.sh
      With NO parallelization extension considered: https://gist.github.com/tkrahn/ef62cfaab678f447ea53ddee09ce0eb2
      Original Extract23 did only for a single 23andMe v3 template.  Here we created an AllSNPs template
      of the merger of all known formats.  Then wrote our own aconv.py module to subset the AllSNPs for each
      desired file format. Returns True if reused existing AllSNPs file; false if created; none if error
    """
    reused = True  # for clarity below (return None on error)

    out_FPB_AllSNPs = f'{wgse.outdir.FPB}_AllSNPs'          # New, v5 true-all form

    # Find needed AllSNPs VCF file to use as a template in variant caller
    refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, type="Microarray")
    if "error" in refVCFtab_qFN or not os.path.isfile(nativeOS(unquote(refVCFtab_qFN))):
        wgse_message("error", "errRefLibMissingFileTitle", True,
                     wgse.lang.i18n["errRefLibMissingAllTable"].replace("{{FILE}}", refVCFtab_qFN))
        return

    # If a valid file exists already; simply return that it is being reused
    if _check_valid_microarray(nativeOS(out_FPB_AllSNPs), unquote(refVCFtab_qFN)):
        return reused

    refgen_qFN = wgse.BAM.Refgenome_qFN
    # missing_refgenome() reports an error if reference genome does not exist; queries user to install if not there
    if wgse.reflib.missing_refgenome(refgen_qFN, parent=parent_window):
        return

    # Setup needed (simplified) file names / paths
    out_FP = f'{wgse.outdir.FP}'
    bamfile = f'"{wgse.BAM.file_FN}"'
    lifthg38 = f'"{wgse.prog_FP}hg38tohg19.py"'
    ploidy_qFN = f'"{wgse.reflib.cma_FP}ploidy.txt"'
    temp_head_qFN = f'"{wgse.reflib.cma_FP}raw_file_templates/head/23andMe_V3.txt"'  # Use 23andMev3 header for AllSNPs

    # Setup needed (simplified) commands (and args to commands)
    # samtools = wgse.samtoolsx_qFN
    bcftools = wgse.bcftoolsx_qFN
    tabix = wgse.tabixx_qFN
    sed = wgse.sedx_qFN
    sort = wgse.sortx_qFN
    cat = wgse.catx_qFN
    zip = wgse.zipx_qFN
    unzip = wgse.unzipx_qFN
    mv = wgse.mvx_qFN
    python = wgse.python3x_qFN
    cpus = wgse.os_threads

    # Setup additional filenames and paths
    temp_called_vcf = f'"{wgse.tempf.FP}AllSNPs_called.vcf.gz"'
    temp_annotated_vcf = f'"{wgse.tempf.FP}AllSNPs_annotated.vcf.gz"'
    temp_result_tab = f'"{wgse.tempf.FP}AllSNPs_result.tab"'
    temp_sorted_result_tab = f'"{wgse.tempf.FP}AllSNPs_result_sorted.tab"'

    # In case DEBUG mode on and Temp directory not cleared out; lets clear these files to make sure not reused
    wgse.tempf.clean_item(temp_called_vcf)        if os.path.exists(temp_called_vcf)        else ""
    wgse.tempf.clean_item(temp_annotated_vcf)     if os.path.exists(temp_annotated_vcf)     else ""
    wgse.tempf.clean_item(temp_result_tab)        if os.path.exists(temp_result_tab)        else ""
    wgse.tempf.clean_item(temp_sorted_result_tab) if os.path.exists(temp_sorted_result_tab) else ""

    commands = ""

    # Rename old vs/v4 CombinedKit to the new CombinedKit_v3 designation; do immediate so can check again in loop
    out_FPB_cmbkitold = f'{wgse.outdir.FPB}_CombinedKit'    # Original name from v3 and v4
    # out_FB_cmbkitold = f'{wgse.BAM.file_FB}_CombinedKit'
    
    cmbkitold_oFN = nativeOS(f'{out_FPB_cmbkitold}.txt')
    cmbkitzip_oFN = nativeOS(f'{out_FPB_cmbkitold}.zip')

    if os.path.exists(cmbkitold_oFN):
        commands += f'{mv} "{out_FPB_cmbkitold}.txt" "{out_FPB_cmbkitold}v3.txt" \n'
    if os.path.exists(cmbkitzip_oFN):
        commands += f'{mv} "{out_FPB_cmbkitold}.zip" "{out_FPB_cmbkitold}v3.zip" \n'
    if os.path.exists(cmbkitzip_oFN) and not os.path.exists(cmbkitold_oFN):
        commands += f'{unzip} -oj "{out_FPB_cmbkitold}v3.zip" -d "{out_FP}" "{out_FPB_cmbkitold}v3.txt" \n'

    # Ready to create AllSNPs command stream

    #  Note: AllSNPs only includes valid SNP calls; not no-call for ones missing or InDels
    #  Note: In v3 added -V indels so we do not generate indel calls and hopefully not indel entries from vendors

    #  Original versus new pileup command; old samtools generates a warning that we should switch to bcftools
    #     {samtools} mpileup -B    -C 50 -r {fchrom} -l {reftab_qFN} -f {refgenome} -hu {bamfile}
    #     {bcftools} mpileup -B -I -C 50 -r {fchrom} -T {reftab_qFN} -f {refgenome} -Ou {bamfile}
    #   Note -B required for Nanopore long read tools; minimal effect on paired-end sequencer output

    # Normal case -- run caller to generate a VCF
    if not parallel:

        commands = (
            f'{bcftools} mpileup -B -I -C 50 -T {refVCFtab_qFN} -f {refgen_qFN} -Ou {bamfile} | '
            f'  {bcftools} call --ploidy-file {ploidy_qFN} -V indels -m -P 0 --threads {cpus} -Oz -o {temp_called_vcf}\n'
        )

    # Parallel DO NOT USE ; INACCURATE (have it on the developers / debug tab)
    # Todo likely could be fixed by adding all the alt-contigs to each list of chromosomes to process.  If know which
    #  alt-contig is tied to each chromosome, then can create proper subsets (which is likely the better idea).
    #  chromosome 1 (and some of the other large ones) could be subsetted to make them smaller for each run.
    else:
        # NOTE: Parallelization does NOT WORK on any platform (tried Win10, Ubuntu, MacOS; Intel, AMD, M1)
        #   Generates about 1% of values in error; especially chr 13-16 and 22 at start. Identical error on all
        #   platforms. Was true for 1K Genome ref models -- bcftools call clearly makes use of alt contigs.
        #   Left it available in the DEBUG tab for further experimentation. Maybe works for HG models?

        commands = (
            f'set +x\n'  # Turn off command echo; otherwise get echo of echo. Comments are not printed unless -v
            f'echo "Generating AllSNPs file from BAM (takes upwards of an hour)" \n'
            f'echo "Parallelizing each chromosome mpileup / variant-call; then will merge the results" \n'
            f'set -x\n'
        )

        # Proper ploidy setting messed up AllSNPs generation for males; so override to avoid warning message
        #  and later error by using special ploidy.txt to give diploid for now.
        # ploidy = "GRCh38" if wgse.BAM.Build == 38 else \
        #          "GRCh37" if wgse.BAM.Build == 37 or wgse.BAM.Build == 19 else "X"
        for chrom in wgse.valid_chromosomes:    # Using internal list; always same chromosome name convention

            chrtmp = f'"{wgse.tempf.FP}{chrom}.vcf.gz"'  # Put each chromosome VCF in temp for later cleanup
            # For GRCh numbering, need to remove "chr" from valid_chromosome name and add T to Mito name
            fchrom = chrom.replace('chr', '').replace("M", "MT") if wgse.BAM.SNTypeC == "Num" else chrom

            commands += (
                f'{bcftools} mpileup -B -I -C 50 -r {fchrom} -T {refVCFtab_qFN} -f {refgen_qFN} -Ou {bamfile} | '
                f'  {bcftools} call --ploidy-file {ploidy_qFN} -V indels -m -P 0 --threads {cpus} -Oz -o {chrtmp} &\n'
            )

        # Wait for all 25 forked processes to finish; then do a subshell cd otherwise shell line exceeds limit
        #  when using 25 fully qualified file names; all files are in temp directory and will get deleted at end
        commands += (
            f'wait\n'
            f'(cd {wgse.tempf.FP}; {bcftools} concat -Oz -o {temp_called_vcf} '
            f' chr[1-9].vcf.gz chr[1-2][0-9].vcf.gz chr[M,X,Y].vcf.gz )\n'
        )

    # Todo more simplication and reduction of intermediate files possible? sed/sort/cat on one line? Avoid more
    #  compression and tabix calls by piping directly?
    commands += (
        f'{tabix} -p vcf {temp_called_vcf}\n'
        f'{bcftools} annotate -Oz -a {refVCFtab_qFN} -c CHROM,POS,ID {temp_called_vcf} > {temp_annotated_vcf}\n'
        f'{tabix} -p vcf {temp_annotated_vcf}\n'
        f'{bcftools} query -f \'%ID\t%CHROM\t%POS[\t%TGT]\n\' {temp_annotated_vcf} -o {temp_result_tab}\n'
        f'{sed} \'s/chr//; s/\tM\t/\tMT\t/g; s/\///; s/\.\.$/--/; s/TA$/AT/; s/TC$/CT/; s/TG$/GT/; s/GA$/AG/; '
        f'    s/GC$/CG/; s/CA$/AC/\' {temp_result_tab} | {sort} -t $\'\t\' -k2,3 -V > {temp_sorted_result_tab}\n'
        f'{cat} {temp_head_qFN} {temp_sorted_result_tab} > "{out_FPB_AllSNPs}.txt"\n'
    )

    # Add liftover command; if needed. Does in place liftover (in the same file)
    if wgse.BAM.Build == 38:
        # refmod = "hg" if wgse.BAM.SNTypeC == "Chr" else "GRCh" if wgse.BAM.SNTypeC == "Num" else "error"
        refmod = "GRCh"     # All microarray files use Numeric naming; so keep the AllSNPs the same
        if refmod == "error":
            wgse_message("error", "errBuildUnkSeqNameTitle", True,
                         wgse.lang.i18n["errBuildUnkSeqNameType"].replace("{{SNType}}", wgse.BAM.SNTypeC))
            return
        commands += f'{python} {lifthg38} "{out_FPB_AllSNPs}" {refmod}\n'

    elif wgse.BAM.Build != 37 and wgse.BAM.Build != 19:
        start = "T2T" if wgse.BAM.Build == 99 else str(wgse.BAM.Build)
        wgse_message("error", "errLiftoverUnsupTitle", True,
                     wgse.lang.i18n["errLiftoverUnsupported"].replace("{{START}}", start).replace("{{END}}", "37"))
        return

    # Final step; create a ZIP file from the TXT file; will check for either on reuse but really want .txt persitent
    # commands += f'{zip} -j "{out_FPB_AllSNPs}.zip" "{out_FPB_AllSNPs}.txt"\n'

    run_bash_script("ButtonAllSNPs", commands, parent=parent_window)

    # Check if now exists; return not reused if so (else None to indicate error and not available)
    if _check_valid_microarray(nativeOS(out_FPB_AllSNPs), unquote(refVCFtab_qFN)):
        return not reused


def cancel_miroarray_window():
    """ Remove Autosomal formats window in preparation for restoring WGSE Main Window """
    from mainwindow import mainwindow_resume   # Have to localize due to import loop

    global microarrayWindow

    try:        # May have been destroyed before this call
        microarrayWindow.destroy()
    except:
        pass
    if __name__ != '__main__':      # for standalone operation
        mainwindow_resume()         # sets up for return to main window


# ***************MAIN PROGRAM (stand-alone)*******************
# Todo standalone program use in development; not been tested of late. Should allow vendors to be specified so no pop-up
if __name__ == '__main__':
    from mainwindow import set_BAM_file, set_outdir       # Need local import to avoid loop

    if len(sys.argv) != 3:
        print("Usage: microarray {bam_file} [outdir]")
        print("       Use fully qualified paths. Outdir is optional. If not specified, the default.")
        print("       is assumed unless specified in the saved settings from a previous run.")
        print("       Command line mode still presents a pop-up selector to specify the formats desired")
        exit()

    """ Microarray start as an independent task (main program) """
    print('Starting WGS Extract Microarray Combined Kit Generator ...')
    wgse.init(interactive=False)

    # Could potentially just get parameters from stored settings JASON file? Or make args optional and call button?
    set_outdir(sys.argv[2], user_set=True)
    if not (wgse.BAM or set_BAM_file(sys.argv[1])):
        exit()      # Already reported issue

    _button_AllSNPs(None)
