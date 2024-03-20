# coding: utf8
#
# Main WIndow / GUI module for the WGS Extract system
# Todo modify from tkinter reliance to another, more modern GUI
# Todo modify to remove button actions so a command line version can be run without the GUI
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    Main window GUI module for WGS Extract.  Everything revolves around the main window and its buttons on
     different panes. Once various support modules are initialized in the main program (wgesextract), the
     main entry point is mainwindow_setup.  Most of the UI and button click action is handled here.
     Notable exception is the microarray generator window pop-up (in microwindow.py)
"""

import os.path      # os.path.isdir, .exists
from functools import partial   # Better than Lambda for some cases
import time         # time.ctime() for results window
import webbrowser   # webbrowser.opennew()
# noinspection PyUnresolvedReferences
import multiqc      # multiqc.run() used only on FastQC output when paired FASTQs to merge into single report

# Original code imported tk after ttk. Thus making tk have priority.
# We know because Tk uses parameter font while Ttk uses styles like CSS; Tk has Message, Ttk only has Label
from tkinter import Tk, Toplevel, filedialog, simpledialog, messagebox, Scrollbar, StringVar    # , Radiobutton
from tkinter import Label, LabelFrame, Text, VERTICAL, HORIZONTAL, LEFT, CENTER, E, W, NW, NSEW, NS, EW, TRUE, SOLID
from tkinter.ttk import Notebook, Frame, Style, Combobox
try:
    from tkmacosx import Button
except ImportError:
    from tkinter import Button

from PIL import ImageTk, Image      # , ImageGrab -- not available in Linux
# noinspection PyUnresolvedReferences,PyPep8Naming
import pyscreenshot as ImageGrab    # Is available in Linux (whereas the original PIL / Pillow is not)

# Local modules from WGSE
from utilities import DEBUG, nativeOS, universalOS, unquote, wgse_message, offset
from utilities import is_legal_path, check_exists, fastq_base
from utilities import is_wgse_outdated, update_release_track, load_wgse_version
from utilities import OutputDirectory, FontTypes

from commandprocessor import run_bash_script, update_release
from bamfiles import BAMFile, BAMContentError, BAMContentErrorFile, BAMContentWarning
from microwindow import button_microarray_window, _button_AllSNPs
from fastqfiles import process_FASTQ

import settings as wgse
font = {}
white = 'snow'          #            : text for low luminance colors (black not visible)
maroon = 'dark violet'  #            : text color, along with red, in normal background
pink = 'pink'           # '#ffb3fe'  : commands less than an hour or two
red = 'hot pink'        # deep pink  : commands more than 2-3 hours
cyan = 'cyan2'          # cyan       : changed from default value
# For text colors use 'red' and '


# Pop-up result windows (main window stored in wgse.window) global to this file
global binstatsWindow, statsWindow, yHaploResWindow, simResWindow, selectLanguageWindow, askRefgenWindow
# global errPopupWindow (utilities), microarrayWindow (microarray)
# global askRefgenWindow (bamfiles), pleaseWaitWindow (commandprocessor)

# Named Buttons (and labels) that are global to this file
global all_action_buttons

global versionLabel, documentationButton, exitButton, titleFileLabel, titleCommandLabel
global outdirLabel, outdirButton, outdirUnselectButton, langSelectLabel, langSelectButton, lreloadButton
global reflibDirectoryLabel, reflibDirectoryButton, tempDirectoryLabel, tempDirectoryButton
global prefServerButton, bamSelectedLabel, bamSelectButton, bamReferenceGenomeLabel
global infoLabelMapAvgReadDepth, infoLabelGender, infoLabelChroms
global bamMapAvgReadDepthLabel, bamGenderLabel, bamChromsLabel, bamFsizeLabel  # bamAverageReadLengthLabel,

global bamIdxstatsButton, bamHeaderButton, bamIndexButton, bamSortButton, bamConvertButton, bamWESButton
global bamRealignButton, bamUnselectButton
global autosomalFormatsButton
global mitoFASTAButton, mitoBAMButton, mitoVCFButton
global yANDmtButton, yOnlyButton, yVCFButton
global haplogroupYButton, haplogroupMtButton, exportUnmappedReadsButton
global bamAlignButton, bamUnalignButton, fastqFastpButton, fastqFastqcButton
global SNPVCFButton, InDelVCFButton, CNVVCFButton, SVVCFButton, AnnotateVCFButton, FilterVCFButton, VarQCButton

# Debug_MODE only buttons, labels
global wslbwaButton, runMicroParallelButton, subsetBAMButton
global maxsettingsLabel, maxmemButton, maxmemLabel, maxthreadButton, maxthreadLabel
global fontsetLabel, fontsizeButton, fontfaceButton, wresetButton, releaseTrackButton, updateProgramButton


def mainwindow_resume():
    """ General function to resume / restore the main window after some user button action(s). """
    update_action_buttons()
    wgse.window.update()
    wgse.window.deiconify()
    wgse.tempf.clean()


def _prepare_to_exit():
    wgse.save_settings()
    if wgse.window:
        wgse.window.withdraw()
    wgse.tempf.list.append(wgse.tempf.oFP)      # Special case; to wipe out the PID created directory
    wgse.tempf.clean(False)                     # Will wipe the temp directory from the list first; so no temp clean


def _button_actual_exit():
    if wgse.window:
        wgse.window.destroy()
    exit()


def button_exit():
    """
    Processing user button (or window close) as main exit out of the program; do any final cleanup.
    Split into two due to Update Release Button to avoid race condition on saving settings
    """
    _prepare_to_exit()
    _button_actual_exit()


def gui_popup(etype, ttitle, tbody):
    """ GUI portion of wgse_message; to keep all tkinter stuff together """

    if etype == "error":
        messagebox.showerror(ttitle, tbody)
    elif etype == "warning":
        messagebox.showwarning(ttitle, tbody)
    elif etype == "info":
        messagebox.showinfo(ttitle, tbody)
    elif etype == "yesno":          # returns True (Yes), False (No)
        return messagebox.askyesno(ttitle, tbody)
    elif etype == "okcancel":       # returns True (OK), False (Cancel)
        return messagebox.askokcancel(ttitle, tbody)
    elif etype == "yesnocancel":    # returns True (Yes), False (No), None (Cancel)
        return messagebox.askyesnocancel(ttitle, tbody)


class ToolTip(object):
    """
    Pop-up / roll-over help messages for buttons, labels and other widgets. Helps alleviate descriptive text in UI.
    """

    def __init__(self, widget):
        self.text = None
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        """Display text in tooltip window"""
        # Todo set geometry?  Or wrapwidth? Or use global text / font setting?
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(TRUE)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


def AddToolTip(widget, text):
    toolTip = ToolTip(widget)

    def enter(event):
        toolTip.showtip(text)

    def leave(event):
        toolTip.hidetip()

    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)


def ask_vcfs_to_process():
    required_ext = (".vcf", ".vcf.gz")
    files = (".vcf", ".gz") if wgse.os_plat == "Darwin" else required_ext   # Continually reintroduced MacOS library bug
    initialdir = wgse.BAM.file_FP if wgse.BAM and wgse.BAM.file_FP else \
                 wgse.outdir.FP if wgse.outdir and wgse.outdir.FP else ""
    vcf_FNs = filedialog.askopenfilenames(parent=wgse.window, initialdir=initialdir,
                                          title=wgse.lang.i18n['SelectVCFFiles'],
                                          filetypes=[(wgse.lang.i18n['FrameVCFFiles'], files), ])

    if not vcf_FNs or len(vcf_FNs) == 0:      # User must have cancelled; simply return as nothing to do
        DEBUG("VCF File Selection Cancelled ...")   # Do nothing
        return ""        # Does not really matter

    DEBUG(f"Selected VCFs (request): {vcf_FNs}")
    VCFs = {}
    for ofile in vcf_FNs:
        if not (ofile and ofile.endswith(required_ext) and os.path.exists(ofile)):
            wgse_message("error", 'VCFFileBadTitle', True, wgse.lang.i18n['VCFFileBad'].replace("{{VCF}}", ofile))
            return ""

        file = universalOS(ofile)
        if not VCFs:    # First one so use it to set a new default output directory if needed
            wgse.outdir.new_default(os.path.dirname(file))

        VCFs[file] = ({'unk'}, set())

    return VCFs


def ask_fastq_to_align():
    required_ext = (".fq", ".fastq", ".fastq.gz", ".fq.gz")
    files = (".fq", ".fastq", ".gz") if wgse.os_plat == "Darwin" else required_ext  # Bug in MacOS ; two dots in extension not allowed
    initialdir = wgse.BAM.file_FP if wgse.BAM    and wgse.BAM.file_FP else \
                 wgse.outdir.FP   if wgse.outdir and wgse.outdir.FP   else ""
    fastq_FNs = filedialog.askopenfilenames(parent=wgse.window, initialdir=initialdir,
                    title=wgse.lang.i18n['SelectFASTQFiles'],
                    filetypes=[(wgse.lang.i18n['FrameFASTQFiles'], files), ])

    if not fastq_FNs or len(fastq_FNs) == 0:      # User must have cancelled; simply return as nothing to do
        DEBUG("FASTQ File Selection Cancelled ...")   # Do nothing
        return ""        # Does not really matter

    for file in fastq_FNs:
        if not (file and file.endswith(required_ext) and os.path.exists(file)):
            wgse_message("error", 'FastqFileBadTitle', True, wgse.lang.i18n['FastqFileBad'].replace("{{FASTQ}}", file))
            return ""

    DEBUG(f"Selected FASTQs (request): {fastq_FNs}")

    # Set default output directory (if not user set yet) to the location the user desires these FASTQs
    wgse.outdir.new_default(os.path.dirname(fastq_FNs[0]))

    return fastq_FNs


def ask_BAM_filename(new=False):
    """
    From Align command; asking for the filename to use in saving the created BAM file. FASTQ files tend to be named
    with complicated forms to simply derive from that. Use saveas instead of open as file does not yet exist.
    """
    title = 'SpecifyBamFile' if new else 'SelectBamFile'
    required_ext = (".bam", ".cram")
    initialdir = wgse.outdir.FP   if wgse.outdir and wgse.outdir.FP   else \
                 wgse.BAM.file_FP if wgse.BAM    and wgse.BAM.file_FP else ""
    file_BS = filedialog.asksaveasfilename(parent=wgse.window, initialdir=initialdir,
                                           title=wgse.lang.i18n[title],
                                           filetypes=[(wgse.lang.i18n['BamFiles'], required_ext), ])

    if not (file_BS and file_BS.endswith(required_ext) and
            (is_legal_path(file_BS) if new else os.path.exists(file_BS))):
        wgse_message("error", 'BAMNameOrExtTitle', False, 'BAMNameOrExt')
        return ""

    if not new:
        wgse.outdir.new_default(os.path.dirname(file_BS))

    return file_BS


global cbRefgenome, askRefgenDone


def ask_refgenome(inBAM=True, window=wgse.window):
    """
    If we cannot figure the reference genome from BAM, ask the user.  Gives them the option to say "unknown".
    This is called from the referencelibrary Class (where the refgenome is determined for a new BAM). But also
    from the (Re)Align button here.  Moved this from referencelibrary class to consolidate GUI code here. As a
    result, we must call back to the referencelibrary to get the list of possible refgenomes to select from
    (as loaded from genomes.csv and possibly trimmed based on known BAM parameters). Unlike other ask_* routines in
    mainwindow, this is a popup selector we generate and not a file selector.  So we need a post pop-up matching
    routine after button selection. We forward the inBAM parameter to the reference library list generation call
    so it can trim the list based on known BAM file parameters that prevent some from being selected.
    Returns a typle of the Refgenome code (string) and the SN count for that reference model (integer)
    """
    global font, askRefgenWindow, cbRefgenome, askRefgenDone
    font = wgse.fonts.table if wgse and wgse.fonts else {'14': ("Times New Roman", 14)}

    # mainWindow.withdraw()

    # Technically, may not be mainwindow that we are forked from; could be anywhere. But have no current window.
    askRefgenWindow = Toplevel(window)
    askRefgenWindow.transient()
    askRefgenWindow.protocol("WM_DELETE_WINDOW", ask_refgenome_cbhandler)   # Trick to clear value
    askRefgenWindow.title(wgse.lang.i18n['SelectReferenceGenome'])
    askRefgenWindow.geometry(offset(window.geometry(), 50))
    askRefgenWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    askRefgenWindow.columnconfigure(0, weight=1)
    askRefgenWindow.rowconfigure(1, weight=1)

    # Setup intro text above list box selector (longer if inBAM and so we could not figure out reference genome)
    textToDisplay = ""
    if inBAM:
        # Failed trying to determine Reference Genome of BAM; so ask user but downselect possible refgenomes
        file = wgse.BAM.file_FB + "_header.txt"
        textToDisplay = wgse.lang.i18n['CouldntDetermineReferenceGenome'].replace("{{Head}}", file) + '\n' + '\n'
        textToDisplay += wgse.lang.i18n['CautionSelectTheCorrectRefGenome'] + '\n' + '\n'

    # Simply ask for a reference genome (e.g. in an Align command) ; if in BAM, may have before and after text
    textToDisplay += wgse.lang.i18n['PleaseSelectReferenceGenome'] + '\n' + '\n'

    if inBAM:
        textToDisplay += '          ' + wgse.lang.i18n['YourBAMParameters'].replace(
            "{{RefGen}}", wgse.BAM.refgenome_str())

    Label(askRefgenWindow, text=textToDisplay, font=font['14'],
          anchor="w", justify="left").grid(row=0, columnspan=2, padx=5, pady=3)

    # Setup listbox dropdown selector
    cbTable = wgse.reflib.available_refgenomes(inBAM)   # table of available ref genomes to select (smaller if inBAM)
    cbRefgenome = StringVar(askRefgenWindow)    # Setup return string variable of drop down table select (cbTable entry)

    askRefgenCB = Combobox(askRefgenWindow, textvariable=cbRefgenome, values=cbTable, state='readonly',
                           height=15, width=70, justify="left")
    askRefgenCB.set('')     # Could default the combobox to a list element e.g. cbTable[0]
    askRefgenCB.bind('<<ComboboxSelected>>', ask_refgenome_cbhandler)
    askRefgenCB.grid(row=1, sticky=NW)

    # We could simply exit the window when the Combobox value is changed; but let's wait for the user to hit a Done
    #  button so they can be assured of what they selcted first
    askRefgenDone = Button(askRefgenWindow, text=wgse.lang.i18n['Done'], font=font['14'], state="disabled",
                           command=ask_refgenome_donehandler)
    askRefgenDone.grid(row=3)

    askRefgenWindow.update()
    askRefgenWindow.grab_set()  # Toplevel equivalent of mainloop
    askRefgenWindow.wait_window(askRefgenWindow)

    # Rertieve user selection and convert back into internal refgenome form (from available selection string)
    new_refgenome = wgse.reflib.key_from_availist(cbRefgenome.get())

    if new_refgenome == "unknown":
        new_refgenome = None

    # If hit "x' to close window without selecting or selected unknown, then new_refgenome is None
    # Otherwise, new_refgenome is the internal key for the selected refgenome string (set by the handler)
    if inBAM:       # Store result in BAM class settings
        wgse.BAM.Refgenome = new_refgenome
        wgse.BAM.RefMito = "Yoruba" if "hg19" in wgse.BAM.Refgenome else \
                           "rCRS" if wgse.BAM.Refgenome else \
                           None
        wgse.BAM.RefgenomeNew = None  # Todo How to set up new RefGenomeNew after the fact? Store in genomes.csv?

    mainwindow_resume()

    return new_refgenome        # Was a function call; so will not return to the main loop


def ask_refgenome_cbhandler(event=None):        # Window exit or combobox selection event (overloaded)
    global askRefgenWindow, cbRefgenome, askRefgenDone

    if event is None:                               # Must be the window close, so clear setting and exit
        cbRefgenome.set("")
        askRefgenWindow.destroy()
    elif cbRefgenome.get() is not None:             # Changed value in combobox so enable Done button
        askRefgenDone.config(state="normal")
    else:                                           # Cleared value from combobox so disable Done button
        askRefgenDone.config(state="disabled")
    # Ignore the value change as we will pick it up when the Done button is hit and we exit the window


def ask_refgenome_donehandler(event=None):      # Done button (only)
    global askRefgenWindow

    askRefgenWindow.destroy()


def update_action_buttons():
    """
    Called to change the status of all main window action buttons; possibly BAM data dependent.

    Newstate is "normal" IF a BAM is set and the output directory is set.  Otherwise, is "disabled"
    Cannot enable Stats nor Header button as it has a Save button which needs output set
    """
    global font

    font = wgse.fonts.table

    # Special buttons that are always active but change to value set; default value initially set
    # Could maybe change out to use StringVar so updated automatically ....
    if wgse.outdir and wgse.outdir.FP and wgse.outdir.user_set:
        outdirButton.configure(text=wgse.outdir.FB)
        outdirButton.configure(font=font['14b'] if wgse.outdir.user_set else font['14'])
        outdirButton.configure(bg=cyan if wgse.outdir.user_set else wgse.defbut_bg)
        outdirUnselectButton.grid()
    else:
        outdirButton.configure(text=wgse.lang.i18n['Default'])
        outdirButton.configure(font=font['14'])
        outdirButton.configure(bg=wgse.defbut_bg)
        outdirUnselectButton.grid_remove()

    reflibDirectoryButton.configure(text=wgse.reflib.FB if wgse.reflib and wgse.reflib.set else
                                         wgse.lang.i18n['Default'])
    tempDirectoryButton.configure(text=wgse.tempf.FB if wgse.tempf and wgse.tempf.set else
                                       wgse.lang.i18n['Default'])
    wslbwaButton.configure(text=wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else
                                wgse.lang.i18n['Inactive'])

    # Relies on global listing "all_action_buttons" that can be set active ("normal") or deactive ("disabled")
    # Remember Stats is automatically run for indexed BAMs
    newstate = "normal" if wgse.outdir and wgse.outdir.FP and wgse.BAM and wgse.BAM.Stats else "disabled"
    for single_button in all_action_buttons:
        single_button.configure(state=newstate)

    # Override BAM File area buttons if BAM and Output Directory are set even if Stats is not yet set
    if wgse.BAM and wgse.outdir and wgse.outdir.FP:
        set_BAM_window_settings()

    # if wgse.outdir and wgse.outdir.FP:                  # Primary input buttons only require outdir be set
    bamSelectButton.configure(state="normal")       # For BAM / CRAM file selection; may be set from saved settings
    bamAlignButton.configure(state="normal")        # For FASTQ; can do Align, fastp, and fastqc as soon as output
    fastqFastpButton.configure(state="normal")      # directory is specified.  It will query for FASTQ files if
    fastqFastqcButton.configure(state="normal")     # not from BAM when BAM specified
    VarQCButton.configure(state="normal")           # Can work off of user specified VCF; not just BAM one
    AnnotateVCFButton.configure(state="normal")     # Can work off of user specified VCF; not just BAM one
    FilterVCFButton.configure(state="normal")       # Can work off of user specified VCF; not just BAM one

    if wgse.os_plat != "Windows":
        wslbwaButton.configure(state="disabled")

    # ******************************************************
    # Todo *** Temporary until we implement the buttons ***
    # SNPVCFButton.grid_remove()
    # InDelVCFButton.grid_remove()
    CNVVCFButton.grid_remove()
    SVVCFButton.grid_remove()
    # VarQCButton.grid_remove()
    # AnnotateVCFButton.grid_remove()
    # FilterVCFButton.grid_remove()
    # ******************************************************

    if newstate == "disabled":
        return
    # else, if newstate == "normal" then Stats have been run and the BAM and Output Directory are both set

    # Special override: Here only if enabled buttons (state = "Normal")
    #   Disable buttons if data is not available (relies on Stats run)
    if wgse.BAM.chrom_types["A"] <= 1:
        autosomalFormatsButton.configure(state="disabled")
        # SNPVCFButton.configure(state="disabled")      # Have Mito and Y VCF buttons; guess could allow for X only?
        # InDelVCFButton.configure(state="disabled")
        # CNVVCFButton.configure(state="disabled")
        # SVVCFButton.configure(state="disabled")
        runMicroParallelButton.configure(state="disabled")

    if wgse.BAM.chrom_types["M"] <= 1:
        mitoFASTAButton.configure(state="disabled")
        mitoBAMButton.configure(state="disabled")
        mitoVCFButton.configure(state="disabled")
        haplogroupMtButton.configure(state="disabled")
        yANDmtButton.configure(state="disabled")

    # If gender is Female, disable Y buttons (note: females have >1 mapped reads on Y!)
    if wgse.BAM.chrom_types["Y"] <= 1 or wgse.BAM.gender == 'Female':
        haplogroupYButton.configure(state="disabled")
        yANDmtButton.configure(state="disabled")
        yOnlyButton.configure(state="disabled")
        yVCFButton.configure(state="disabled")

    if wgse.BAM.chrom_types["*"] <= 1:
        exportUnmappedReadsButton.configure(state="disabled")

    if wgse.BAM.Build in [36, 99]:      # Build 36 and 99 missing support files (yVCFButton simply report error)
        autosomalFormatsButton.configure(state="disabled")
        haplogroupYButton.configure(state="disabled")
        haplogroupMtButton.configure(state="disabled")
        # SNPVCFButton.configure(state="disabled")
        runMicroParallelButton.configure(state="disabled")

    if wgse.BAM.RefMito != "rCRS":
        haplogroupMtButton.configure(state="disabled")

    # SNPVCFButton.configure(state="disabled")
    # InDelVCFButton.configure(state="disabled")
    CNVVCFButton.configure(state="disabled")
    SVVCFButton.configure(state="disabled")


def clear_BAM_window_settings():
    """ Called after clearing a BAM to clear all the window labels of wgse.BAM.* settings. """

    if not wgse.window:
        return

    # Title bar filename label
    titleFileLabel.configure(text="")

    # Reset button to initial Select message
    bamSelectButton.configure(text=wgse.lang.i18n['SelectBamFile'])
    # if not (wgse.outdir and wgse.outdir.FP):          # Setup default outdir for BAM now; so do not disable
    #    bamSelectButton.configure(state="disabled")

    bamReferenceGenomeLabel.configure(text="")

    # Unaligned BAMs will change the left label; so reset these to default for Aligned BAMs
    infoLabelMapAvgReadDepth.configure(text=wgse.lang.i18n['MapAvgReadDepthNoN'] + ":")
    bamMapAvgReadDepthLabel.configure(text="")

    infoLabelGender.configure(text=wgse.lang.i18n['Gender'] + ":")
    bamGenderLabel.configure(text="")

    infoLabelChroms.configure(text=wgse.lang.i18n['Chroms'] + ":")
    bamChromsLabel.configure(text="")

    bamFsizeLabel.configure(text="")

    bamIdxstatsButton.configure(state="disabled")
    bamHeaderButton.configure(state="disabled")
    bamSortButton.grid_remove()
    bamIndexButton.grid_remove()
    bamConvertButton.grid_remove()
    bamWESButton.grid_remove()
    bamRealignButton.grid_remove()
    bamUnselectButton.grid_remove()


def set_BAM_window_settings():
    """
    Called after processing a new BAM to (re)set window text labels dependent on wgse.BAM.* settings
     Also called when simply index / unindex a BAM as it changes the Indexed state.
     Note that sort / unsort, to CRAM / BAM, and realign buttons all cause a new BAM to be set which then call this.
     Stats button may cause delayed idxstats to run (for CRAM, unindexed BAM) so button_close also calls here
       button_close reused by other results windows so causes a few extra calls but no harm
    """

    if not wgse.BAM:
        clear_BAM_window_settings()     # Just to make sure; but likely OK to just return
        return

    if not wgse.window:
        return

    # For title bar filename (so visible no matter what the tab setting; as a reminder)
    titleFileLabel.configure(text=wgse.BAM.disp_FBS)

    # Changed from FB to FBS to show CRAM file extension; use size limited one to avoid overruns in display
    bamSelectButton.configure(text=wgse.BAM.disp_FBS)
    # if wgse.outdir and wgse.outdir.FP:        # Always display BAM select button now
    bamSelectButton.configure(state="normal")
    bamReferenceGenomeLabel.configure(text=wgse.BAM.refgenome_str())

    if wgse.BAM.Stats:      # Values to report as "blank" when no Stats run to fill them in yet
        bamMapAvgReadDepthLabel.configure(text=f'{wgse.BAM.mapped_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}')
        bamGenderLabel.configure(
                text=wgse.lang.i18n.get(wgse.BAM.gender, wgse.BAM.gender))  # use value if key error
    elif wgse.BAM.Sorted is None:   # Unaligned BAM; report different information from normal Stats page
        read_type = wgse.lang.i18n['PairedEnd'] if wgse.BAM.ReadType == "Paired" else \
            wgse.lang.i18n['SingleEnd'] if wgse.BAM.ReadType == "Single" else ""
        read_length_str = (f'{wgse.BAM.avg_read_length:,.0f} {wgse.lang.i18n["bp"]}, '
                           f'{wgse.BAM.avg_read_stddev:,.0f} {wgse.lang.i18n["StdDev"]}, {read_type}')
        infoLabelMapAvgReadDepth.configure(text=wgse.lang.i18n['AverageReadLength'] + ":")
        bamMapAvgReadDepthLabel.configure(text=read_length_str)

        infoLabelGender.configure(text=wgse.lang.i18n['Sequencer'] + ":")
        bamGenderLabel.configure(text=wgse.BAM.Sequencer if wgse.BAM.Sequencer else "Unknown")

    else:   # No BAM information to report
        bamMapAvgReadDepthLabel.configure(text="")
        bamGenderLabel.configure(text="")

    bamChromsLabel.configure(text=wgse.BAM.chrom_types_str())   # Now called File Content: (is blank if no info)
    bamFsizeLabel.configure(text=wgse.BAM.filestatus_str())     # Now called File Stats: as has sort and index status

    stats_run = wgse.BAM.Stats or (wgse.BAM.Indexed and wgse.BAM.file_type == "BAM")
    bamIdxstatsButton.configure(bg=wgse.defbut_bg if stats_run else pink)
    bamIdxstatsButton.configure(state="disabled" if wgse.BAM.Sorted is None else "normal")
    bamIdxstatsButton.grid()     # Always show
    bamHeaderButton.configure(state="normal")
    bamHeaderButton.grid()       # Always show

    # For Sort and Index buttons that serve as status indicators also; disable if already done and not DEBUG mode
    #  otherwise, if DEBUG mode, always enabled but a toggle button
    if wgse.DEBUG_MODE and wgse.BAM.Sorted:   # For developers only
        bamSortButton.configure(state="normal", text=wgse.lang.i18n["Unsort"], bg=wgse.defbut_bg, command=_button_unsort_BAM)
    elif wgse.BAM.Sorted is None:
        bamSortButton.configure(state="disabled", text=wgse.lang.i18n["Unaligned"], bg=wgse.defbut_bg, command=_button_unsort_BAM)
    elif wgse.BAM.Sorted:
        bamSortButton.configure(state="disabled", text=wgse.lang.i18n["Sorted"], bg=wgse.defbut_bg, command=button_sort_BAM)
    else:   # BAM is aligned and not coordinate sorted; so display the Sort button
        bamSortButton.configure(state="normal", text=wgse.lang.i18n["Sort"], bg=pink, command=button_sort_BAM)
    bamSortButton.grid()        # Make button visible because BAM File is set now

    if wgse.DEBUG_MODE and wgse.BAM.Indexed:   # For developers only
        bamIndexButton.configure(state="normal", text=wgse.lang.i18n["Unindex"], bg=wgse.defbut_bg, command=_button_unindex_BAM)
    elif wgse.BAM.Indexed:
        bamIndexButton.configure(state="disabled", text=wgse.lang.i18n["Indexed"], bg=wgse.defbut_bg, command=button_index_BAM)
    else:   # BAM is not indexed
        bamIndexButton.configure(state="normal", text=wgse.lang.i18n["Index"], bg=pink, command=button_index_BAM)
    bamIndexButton.grid()       # Make button visible because BAM File is set now

    # Additional buttons that are initially invisible and should now be made visible
    if wgse.BAM.file_type == "BAM":
        bamConvertButton.configure(state="normal", text=wgse.lang.i18n["ToCRAM"], command=button_BAM_to_CRAM)
    elif wgse.BAM.file_type == "CRAM":
        bamConvertButton.configure(state="normal", text=wgse.lang.i18n["ToBAM"], command=button_CRAM_to_BAM)
    bamConvertButton.grid()  # Make button visible

    bamWESButton.configure(state="normal", text=wgse.lang.i18n['ToPozBAM' if wgse.BAM.Yonly else 'ToWESBAM'])
    bamWESButton.grid()

    bamRealignButton.configure(state="normal")
    bamRealignButton.grid()     # Make button visible

    bamUnselectButton.configure(state="normal")
    bamUnselectButton.grid()     # Make button visible


def button_select_BAM_file():
    """
        Processing user button to select a (new) BAM File. Call internal set_BAM_file with new file selected.
    """
    # Ask user for BAM file to load; must be mainWindow button click to set BAM
    initialdir = nativeOS(wgse.BAM.file_FP) if wgse.BAM and wgse.BAM.file_FP else \
                 wgse.outdir.oFP if wgse.outdir and wgse.outdir.oFP else ""
    BAM_FN = filedialog.askopenfilename(parent=wgse.window, initialdir=initialdir,
                title=wgse.lang.i18n['SelectBamFile'],
                filetypes=((wgse.lang.i18n['BamFiles'], (".bam", ".cram", ".sam")),))
    if BAM_FN:  # User selected new file name; so process it
        set_BAM_file(BAM_FN)
    else:       # User must have cancelled; simply return as nothing to do
        DEBUG("BAM/CRAM/SAM File Selection Cancelled ...")   # Do nothing; do NOT unselect current one
    mainwindow_resume()


def button_unselect_BAM_file():
    if wgse.BAM:
        wgse.BAM = None
        if not wgse.outdir.user_set:
            button_unselect_outdir()
        set_BAM_window_settings()           # Sets window labels of BAM file settings
        wgse.save_settings()                # Changed BAM file; so save settings
    mainwindow_resume()


def set_BAM_file(BAM_FN):
    """
        Initiates basic processing to (change/set) BAM file.  Special handling as can have errors deep in the
        call stack. So use  Try...Except with custom error class here. Handle and catch all known errors at this level.
        Otherwise, no general "except:"-ion to catch all errors (simply pass on coding errors).
        This can also be called internally to replace (or set new) a BAM file.  So must diligently handle
        checking and setup before doing actual replacement.
    """
    # Todo change to wgse.BAM.change(); but what if no BAM ever set yet? Do init first then call change?
    DEBUG(f"Set BAM/CRAM/SAM (request): {BAM_FN}")

    # Empty BAM file name; error message and return as nothing to do
    if not BAM_FN:
        wgse_message("error", 'InvalidBAMTitle', True,
                     wgse.lang.i18n['errBAMFile'].replace("{{BAMF}}", "none") + "(blank)")
        return False

    # Shortcut if not changing the currently set BAM file
    if wgse.BAM and wgse.BAM.file_FN == BAM_FN:
        return True

    # Prepare to restore existing BAM (and default outdir) if new one fails; even if BAM not yet set and entry is None
    save_orig_BAM = wgse.BAM
    save_orig_outdir = wgse.outdir

    wgse.BAM = None
    if wgse.outdir and wgse.outdir.FP and not wgse.outdir.user_set:
        wgse.outdir = OutputDirectory()     # Setup new, blank class instance IF no current outdir or not user_set

    try:
        wgse.BAM = BAMFile(BAM_FN)      # Process the BAM by opening file, reading header and setting up Class instance
    except BAMContentErrorFile as err:
        # Identical content at this point to normal error as file was added into reason already
        if err.reason != "RefGen":      # Special to indicate error already reported
            wgse_message("error", 'InvalidBAMTitle', True,
                         wgse.lang.i18n['errBAMFile'].replace("{{BAMF}}", BAM_FN) + err.reason)
    except BAMContentError as err:      # Error setting up new BAM; restore old one and report issue
        wgse_message("error", 'InvalidBAMTitle', True,
                     wgse.lang.i18n['errBAMFile'].replace("{{BAMF}}", BAM_FN) + err.reason)
    except BAMContentWarning as warn:   # Warning on setting up new BAM; simply issue message
        wgse_message("warning", 'FrameBAM', True,
                     wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) + warn.reason)

    else:                               # BAM was setup; check for other issues to warn about
        if not (wgse.BAM.Sorted and wgse.BAM.Indexed) and wgse.BAM.Sorted is not None:
            # Warn about not collecting stats on initial BAM load if not sorted and/or indexed (ignore if Unaligned)
            wgse_message("warning", 'FrameBAM', True,
                         wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) +
                         wgse.lang.i18n['warnBAMNoStatsNoIndex'])
        elif wgse.BAM.file_type == "CRAM" and not wgse.BAM.Stats:  # Do not need to warn about both; one is enough
            # Warn about need to collect stats on CRAM file; if we did not pick up stats from a previous run
            wgse_message("warning", 'FrameBAM', True,
                         wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) +
                         wgse.lang.i18n['warnCRAMNoStats'])

    # If still not set, restore old settings (BAM may be None)
    if not wgse.BAM:
        wgse.BAM = save_orig_BAM            # Replaces with previous class instances (garbage collect old ones)
        wgse.outdir = save_orig_outdir

    # BAM file did not change (exited earlier if it was the same so how did we get here?)
    if wgse.BAM and wgse.BAM.file_FN != BAM_FN:
        return False

    set_BAM_window_settings()           # Sets window labels of BAM file settings
    wgse.save_settings()                # Changed BAM file; so save settings
    return True


def button_stats_BAM():
    """
    Processing user button to show detailed BAM Stats. We have likely already run samtools_stats.
    So maybe just displaying the stored table.
    """

    # There are three stats runs.  Each possibly 30-60 minutes long.  We delay long stats runs until the user
    #  requests it directly with a button click.  But, if they were run previously and we find the file is still
    #  available, we will go ahread and process the file now.  Each stats call determines if it was run before
    #  and/or can run quickly; unless being forced by a button click directly.

    # Stats Button only enabled if wgse.BAM exists; but Stats may not have been calculated yet so force it now ...
    try:
        wgse.BAM.get_samtools_idxstats(button_directly=True)
    except Exception as err:      # Error generating idxstats file
        wgse_message("error", 'errNoFileExists', True, err.reason)
        return

    if wgse.BAM.Stats is False:
        DEBUG("*** Failed to create the idxstats file.")
        return  # This should never have occurred

    # So do not force now (as if a button hit directly). Instead, determine if there is a file from a previous run and
    #  load now if so.
    if wgse.BAM.Primary:   # Unmapped or Alt Contig only BAMs cannot have coverage
        wgse.BAM.get_bincvg_stats(wgse.window, "WGS", button_directly=False)    # button=false only reads stats if found
        wgse.BAM.get_bincvg_stats(wgse.window, "WES", button_directly=False)    # - does not create if not found

    # Create the pop-up window to display the stats
    result_window_stats()


def result_window_stats():
    """
    Create and setup main Stats window. Wait until closed.
    This is the only window (beside the main) that has buttons for creating sub-windows.
    """
    global statsWindow, font

    covstats = wgse.BAM.coverage
    WESstats = wgse.BAM.coverage_WES

    font = wgse.fonts.table
    fontti = font['18b']        # Title font for each frame
    fonth1 = font['14b']        # Header font
    fontb1 = font['14']         # Body font
    fontb2 = font['13']         # Special for by-chr body only
    fontb3 = font['12']         # special for date/time/prog_version; by name table entries

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    statsWindow = Toplevel(wgse.window)
    statsWindow.transient(wgse.window)
    statsWindow.title(wgse.lang.i18n['BamFileStatistics'])
    statsWindow.geometry(wgse.bamstats_winsize + offset(wgse.window.geometry(), 50))
    statsWindow.protocol("WM_DELETE_WINDOW", lambda win=statsWindow: button_close(win))
    statsWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    statsWindow.columnconfigure(0, weight=1)
    statsWindow.rowconfigure(0, weight=1)

    # Setup by-Chromosome Frame (main table, left)
    byChromFrame = LabelFrame(statsWindow, text=wgse.lang.i18n['FrameByRefSeqName'], font=fontti)
    byChromFrame.grid(row=0, column=0, padx=5, pady=2)
    byChromFrame.columnconfigure(0, weight=1)
    byChromFrame.rowconfigure(0, weight=3)

    # Setup Summary Frame (column, right; with Close button)
    summaryFrame = LabelFrame(statsWindow, text=wgse.BAM.disp_FBS, font=fontti)
    summaryFrame.grid(row=0, column=1, padx=5, pady=2)

    # ---------------------------------------------------------------------------------------------------------------
    # by Reference Sequence summary table at left

    # Header for by-Chromosome table
    table_header1 = ['RefSeqName', 'LenInModel', 'NCntInModel', 'ReadSegsMap', 'MappedGbases',
                     'MapAvgReadDepthNoNShort', 'Coverage']  # Moved Coverage button to summary area; label still here
    for col in range(len(table_header1)):   # Double row in grid using rowspan;  2x high using \n in label text
        Label(byChromFrame, text=wgse.lang.i18n[table_header1[col]],
              font=fonth1).grid(column=col, row=0, rowspan=2, padx=5, pady=0)
    # Button for Coverage removed from column heading in v4.
    # Changed WEScvg and Coverage to Bin Coverage commands with separate pop-up by chromosome.

    # Body for by-Chromosome table; including total row with possible Other row trailing after
    for row in range(len(wgse.BAM.stats_chroms)):
        if wgse.BAM.stats_chroms[row][0] == 'G':        # Do not display Grand Total row
            continue
        total_row = wgse.BAM.stats_chroms[row][0] == 'T'
        asterisk = "*" if "Y" in wgse.BAM.stats_chroms[row][0] or \
                          ("X" in wgse.BAM.stats_chroms[row][1] and wgse.BAM.gender == "Male") else None
        #          '+' if "O" in wgse.BAM.stats_chroms[row][0]
        width = 0;  val = ""
        # First column is key with real start in 2nd; but want to add blank for missing last Coverage column
        for col in range(len(wgse.BAM.stats_chroms[row])-1):
            index = col+1
            if col == 0:                            # Chromosome / Sequence label
                val = '{:>5}'.format(wgse.BAM.stats_chroms[row][index])
                width = 5
            if col == 1:                            # Use locale specific grouping (matters if K OR actual; not M)
                large = int(wgse.BAM.stats_chroms[row][index]) > 10 ** 6  # When unmapped is too large for unloc segs
                val = '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                       wgse.lang.i18n["Thous"]) if not total_row and not large else \
                      '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -6)/(10.0 ** 6),
                                        wgse.lang.i18n["Mill"])
                width = 7                           # Room for appended " K"; possible billion number in Total
            elif col in [2, 3]:                     # Use locale specific grouping (matters if K OR actual; not M)
                large = int(wgse.BAM.stats_chroms[row][index]) > 10 ** 6    # When unmapped is too large for unloc segs
                val = '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                       wgse.lang.i18n["Thous"]) if not total_row and not large else \
                      '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -6)/(10.0 ** 6),
                                        wgse.lang.i18n["Mill"])
                width = 7                           # Room for appended " K"; common hundred million max value
            elif col == 4:                          # Process float gbases; skip '*'
                # '*' line may have null strings? (historically)
                val = '{:6.2f}'.format(float(wgse.BAM.stats_chroms[row][index])) \
                    if wgse.BAM.stats_chroms[row][index] else '0.00'
                width = 5
            elif col == 5:                          # Avg Read Depth (add 'x')
                xstar = f'{wgse.lang.i18n["xTimes"]}' + (asterisk if asterisk else "")
                large = int(wgse.BAM.stats_chroms[row][index]) > 1000       # When Mito is too large for Mapped ARDepth
                val = '{:5n}{} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                          wgse.lang.i18n["Thous"], xstar) if large else \
                      '{:6n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), 0), xstar)
                width = 5
            elif col == 6:                          # Breadth of Coverage (extra run / button)
                val = '{:6.2f} %'.format(round(wgse.BAM.stats_chroms[row][index] * 100, 2)) if covstats else ""
                width = 8       # mtDNA is always 100% so need the room
            Label(byChromFrame, text=val, font=fonth1 if total_row else fontb2,
                  width=width, anchor="e").grid(column=col, row=row+2, padx=0, pady=0)      # Was 10, 0

    # ---------------------------------------------------------------------------------------------------------------
    # Overall stats summary Frame at right

    # Summary Frame top section: Mapped / Raw differentiated stats
    crow = 0
    Label(summaryFrame, text=wgse.lang.i18n['Mapped'].upper(), font=fonth1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=wgse.lang.i18n['Raw'].upper(), font=fonth1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageReadDepthNoN'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{wgse.BAM.mapped_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}',
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    labeltext = 'AverageReadDepthPoz' if wgse.BAM.Yonly else 'AverageReadDepthWES'
    Label(summaryFrame, text=wgse.lang.i18n[labeltext], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    if WESstats:            # Summary row WSE stats ready
        Label(summaryFrame, text=f'{wgse.BAM.mapped_avg_read_depth_WES} {wgse.lang.i18n["xTimes"]}',
              font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
        Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_depth_WES} {wgse.lang.i18n["xTimes"]}',
              font=fontb1).grid(row=crow, column=3, padx=0, pady=0)
    crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Gbases'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=str(wgse.BAM.mapped_gbases), font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=str(round(wgse.BAM.raw_gigabases, 2)),
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['TotalSegs'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.mapped_segs_read / (10.0 ** 6), 0))} {wgse.lang.i18n["Mill"]}',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.raw_segs_read / (10.0 ** 6), 0))} {wgse.lang.i18n["Mill"]}',
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['NumReads'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.mapped_reads_percent * 100, 0))} %',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text="100 %",
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text="", font=fontb1).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    # Summary Frame middle section: Coverage buttons / stats
    if wgse.BAM.Primary:   # Only display coverage button and stats if primaries exist (could be unmapped only BAM)
        bg = pink if not covstats else wgse.defbut_bg
        Button(summaryFrame, text=wgse.lang.i18n['CoverageShort'], justify=LEFT, bg=bg,
               command=lambda win=statsWindow: button_statsbin_BAM(win, "WGS"),
               font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
        if covstats:
            Label(summaryFrame, text='{:8.4f} %'.format(round(wgse.BAM.coverage * 100, 4)),
                  font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0)
        crow += 1

        bg = pink if not WESstats else wgse.defbut_bg
        labeltext = 'CoveragePozShort' if wgse.BAM.Yonly else 'CoverageWESShort'
        Button(summaryFrame, text=wgse.lang.i18n[labeltext], justify=LEFT, bg=bg,
               command=lambda win=statsWindow: button_statsbin_BAM(win, "WES"),
               font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
        if WESstats:
            Label(summaryFrame, text='{:8.4f} %'.format(round(wgse.BAM.coverage_WES * 100, 4)),
                  font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1
        crow += 1

    # Summary Frame bottom section: Misc values

    # Prepare some strings first
    read_type = wgse.lang.i18n['PairedEnd'] if wgse.BAM.ReadType == "Paired" else \
        wgse.lang.i18n['SingleEnd'] if wgse.BAM.ReadType == "Single" else ""

    read_length_str = (f'{wgse.BAM.avg_read_length:,.0f} {wgse.lang.i18n["bp"]} {chr(0xB1)} '
                       f'{wgse.BAM.avg_read_stddev:,.0f} , {read_type}')

    if wgse.BAM.ReadType == "Paired":
        insert_size_str = (f'{wgse.BAM.insert_size:,.0f} {wgse.lang.i18n["bp"]} {chr(0xB1)} '
                           f'{wgse.BAM.insert_stddev:,.0f}')
    else:
        insert_size_str = f'({wgse.lang.i18n["SingleEnd"]})'

    map_quality_str = f'Q{wgse.BAM.mapq:.0f} {chr(0xB1)} {wgse.BAM.mapq_stddev:,.0f}'
    base_quality_str = f'{wgse.BAM.pup_q30:.0%} >Q30 ; {wgse.BAM.pup_q20:.0%} >Q20'


    Label(summaryFrame, text=wgse.lang.i18n['RefModel'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.refgenome_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageReadLength'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=read_length_str,
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageInsertSize'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=insert_size_str,
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['DupsPercent'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{wgse.BAM.dups_percent:.1%}',
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['QualityMap'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=map_quality_str,
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['QualityBase'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=base_quality_str,
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Chroms'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.chrom_types_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Gender'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.lang.i18n[wgse.BAM.gender] if wgse.BAM.gender else "",
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Sequencer'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.Sequencer if wgse.BAM.Sequencer else "Unknown",
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['FileStats'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.filestatus_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Button(summaryFrame, text=wgse.lang.i18n['SaveWindow'], font=fonth1,
           command=lambda win=statsWindow: button_save(win, "stats")).grid(row=crow, column=0, padx=5, pady=2)
    Button(summaryFrame, text=wgse.lang.i18n['CloseWindow'], font=fonth1,
           command=lambda win=statsWindow: button_close(win)
           ).grid(row=crow, column=1, columnspan=3, padx=5, pady=2); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['StatsNote'], font=fontb3
          ).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=fontb3).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    # Display and wait for Close button
    statsWindow.update()
    statsWindow.grab_set()          # Toplevel equivalent of mainloop
    statsWindow.wait_window()


def result_window_binstats(window, bamtype="WGS"):
    """
    A subwindow to the Stats window for displaying the detailed Bin Coverage stats. Either WGS or WES table by sequence.
    Note that WES uses the TotalBC sequence length and not the reference model value (which includes Ns).
    Similar to result_window_stats() call.
    """
    global binstatsWindow, font

    if (bamtype == "WGS" and not wgse.BAM.coverage) or \
       (bamtype == "WES" and not wgse.BAM.coverage_WES):
        # Todo generate error pop-up instead of silent return?  Call to display stats that are not available.
        return

    font = wgse.fonts.table
    fontti = font['18b']  # Title font for each frame
    fonth1 = font['14b']  # Header font
    fontb1 = font['14']  # Body font
    fontb2 = font['13']  # Special for by-chr body only
    fontb3 = font['12']  # special for date/time/prog_version; by name table entries

    stats_bin = wgse.BAM.stats_bin if bamtype == "WGS" else wgse.BAM.stats_binwes

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    binstatsWindow = Toplevel(window)
    binstatsWindow.transient(window)
    binstatsWindow.title(wgse.lang.i18n['BamFileStatistics'])
    binstatsWindow.geometry(wgse.bamstats_winsize + offset(window.geometry(), 50))
    binstatsWindow.protocol("WM_DELETE_WINDOW", lambda win=binstatsWindow: button_close(win))
    binstatsWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    binstatsWindow.columnconfigure(0, weight=1)
    binstatsWindow.rowconfigure(0, weight=1)

    # Setup by-Chromosome Frame (main table, left)
    frame_name = 'FrameBinCov' if bamtype == "WGS" else 'FrameBinCovWES'    # Todo account for POZ vs WES
    byChromFrame = LabelFrame(binstatsWindow, text=wgse.lang.i18n[frame_name], font=fontti)
    byChromFrame.grid(row=0, column=0, padx=5, pady=2)
    byChromFrame.columnconfigure(0, weight=1)
    byChromFrame.rowconfigure(0, weight=3)

    # Summary frame is just buttons here (right)
    summaryFrame = LabelFrame(binstatsWindow, text=wgse.BAM.disp_FBS, font=fontti)
    summaryFrame.grid(row=0, column=1, padx=5, pady=2)

    # Header for by-Chromosome table (removed 'bin0' column as superfluous)
    table_header1 = ['RefSeqName', 'LenInModel', 'bin1-3', 'bin4-7', 'bin8-', 'bin1-']
    for col in range(len(table_header1)):
        Label(byChromFrame, text=wgse.lang.i18n[table_header1[col]],
              font=fonth1).grid(column=col, row=0, rowspan=2, padx=5, pady=0)  # Single row but 2x high using \n in text

    # Body for by-Chromosome table; including final total row.  The stats_bin[wes] table columns:
    # chromnum, chromosome, total_bases_cnt, ncount, zero_bases_cnt, bin1_3_cnt, bin4_7_cnt, bin8_up_cnt, nz_bases_cnt
    for row in range(len(stats_bin)):
        if stats_bin[row][0] == 'G':        # Do not display Grand Total row
            continue
        total_row = stats_bin[row][0] == 'T'
        asterisk = "*" if "Y" in wgse.BAM.stats_chroms[row][0] and wgse.BAM.gender == "Female" else None
        #          "+" if "O" in wgse.BAM.stats_chroms[row][0] else None

        # First column is key with real name in 2nd; but will use key (generic name) for consistency
        for col in range(1, len(stats_bin[row])):       # Start with col 1 and not 0
            if col == 1:  # Chromosome / Sequence label
                val = '{:>5}'.format(stats_bin[row][1])
                width = 5
            elif col == 2:  # Use locale specific grouping (matters if K OR actual; not M)
                total_bases_cnt = int(stats_bin[row][2] - stats_bin[row][3])
                large = total_bases_cnt > 10 ** 6  # When unmapped is too large for unloc segs
                val = '{:5n} {}'.format(round(total_bases_cnt, -3) / 1000,
                            wgse.lang.i18n["Thous"]) if not total_row and not large else \
                      '{:5n} {}'.format(round(total_bases_cnt, -6) / (10.0 ** 6),
                            wgse.lang.i18n["Mill"])
                width = 7  # Room for appended " K"; possible billion number in Total
            elif col in [3, 4]:  # Skip ncount and zero count columns
                continue
            else:   # Must be [5, 6, 7, 8] All percentage (fractional) values
                val = '{:6.2f} %'.format(round(stats_bin[row][col] * 100, 2))
                val += asterisk if asterisk and col == 8 else ""
                width = 8   # Might be 100%
                col -= 2    # Reflect skipping column 3 & 4 in original table -- N value count)

            Label(byChromFrame, text=val, font=fonth1 if total_row else fontb2,
                  width=width, anchor="e").grid(column=col-1, row=row + 2, padx=0, pady=0)  # Was 10, 0

    Label(summaryFrame, text=wgse.lang.i18n['Gender'], justify=LEFT,
          font=fonth1).grid(row=0, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.lang.i18n[wgse.BAM.gender] if wgse.BAM.gender else "",
          font=fontb1).grid(row=0, column=1, padx=0, pady=0)

    wtype = "binstats" if bamtype == "WGS" else "wesstats"      # Todo account for POZ instead of WES in Y only
    Button(summaryFrame, text=wgse.lang.i18n['SaveWindow'], font=fonth1,
           command=lambda win=binstatsWindow, wt=wtype: button_save(win, wt)).grid(row=1, column=0, padx=5, pady=2)
    Button(summaryFrame, text=wgse.lang.i18n['CloseWindow'], font=fonth1,
           command=lambda win=binstatsWindow: button_close(win)).grid(row=1, column=1, padx=5, pady=2)

    Label(summaryFrame, text=wgse.lang.i18n['StatsNote'],
          font=fontb3).grid(row=2, column=0, columnspan=2, padx=0, pady=0)

    Label(summaryFrame, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=fontb3).grid(row=3, column=0, columnspan=2, padx=0, pady=0)

    # Display and wait for Close button
    binstatsWindow.update()
    binstatsWindow.grab_set()  # Toplevel equivalent of mainloop
    binstatsWindow.wait_window()


'''
def button_statsCov_BAM(window):    # DEPRECATED: See button_statsbin_BAM() called from  result_window_stats() directly
    """
    A separate button just to get the coverage per SN.  Takes another 30 minutes so let user choose to add the column.
    Initiated by hitting button in Coverage column of main stats page.  Will replace that page with a new one once read.
    IMPORTANT: Takes the covbases column 5 only.  Uses the stored model LN and lookup N to get a better model length
    to determine the actual coverage from.  Makes a big difference on the Y chromosome!
    """
    wgse.BAM.get_coverage_stats(window, button_directly=True)   # Call wgse.BAM routine to run samtools coverage command

    try:
        window.destroy()                # Remove the current stats window (left up while calculating the new one)
    except:
        pass
    button_stats_BAM()     # Call the button to regenerate the stats window but with the Coverage data now available


def button_statsWSE_BAM(window):    # DEPRECATED: See button_statsbin_BAM called from result_window_stats() directly
    """
    A separate button just to calculate the WES average read depth and coverage of sample. Takes another 30 minutes
    so let user choose when they want to add data.  Initiated by hitting the Calculate button next to the WES
    label / title.  Will replace current stats page with new one once done and fill in set values
    """
    wgse.BAM.get_WEScvg_stats(window, button_directly=True)       # Call wgse.BAM routine to run samtools depth on WES region

    try:
        window.destroy()                # Remove the current stats window (left up while calculating the new one)
    except:
        pass
    button_stats_BAM()     # Call the button to regenerate the stats window but with WSE data now available
'''


def button_statsbin_BAM(window, bamtype="WGS"):    # button local to result_window_stats() window
    """
    Generate BinCov stats on primaries and Mito using samtools depth (if not done so already) and then display
    Use WES BED file if scantype set to WES
    """
    stats_found = (bamtype == "WGS" and wgse.BAM.coverage) or (bamtype == "WES" and wgse.BAM.coverage_WES)
    if not stats_found:
        wgse.BAM.get_bincvg_stats(window, bamtype, button_directly=True)

    result_window_binstats(window, bamtype)

    if not stats_found:        # New stats will update the main Stats Window so regenerate
        try:
            window.destroy()   # Remove the current stats window (left open while calculating the new one)
        except:
            pass
        button_stats_BAM()     # Call the button to (re)generate the stats window with coverage data now available


def button_save(window, wtype):
    """ Processing user button to save BAM Stats windows """
    # Find Upper Left and Lower Right corners of Stats Window
    ulx = window.winfo_rootx();  uly = window.winfo_rooty()
    lrx = window.winfo_width();  lry = window.winfo_height()
    ImageGrab.grab().crop((ulx, uly, ulx+lrx, uly+lry)).convert("RGB").save(f'{wgse.outdir.oFPB}_{wtype}.jpg')


def button_close(window, top=True):
    """ Processing user button to close sub-window like BAM stats or some result """
    try:
        window.destroy()
    except:
        pass
    if top:
        set_BAM_window_settings()  # May have updated BAM stats by clicking the stats button, so ...
        mainwindow_resume()


def button_header_BAM():
    """ Display BAM/CRAM/SAM header in a pop-up results window and write text file """

    header_oFN = f'{wgse.outdir.oFP}{wgse.BAM.file_FB}_header.txt'
    with open(header_oFN, "wb") as f:          # 15 Mar 2020 REH Changed to binary write to keep Unix \n format
        f.write(wgse.BAM.Header.encode())

    result_window_simple(wgse.BAM.file_FBS + " Header", wgse.BAM.Header, "header", '12', top=True, wwrap="none")


def _adjust_mem_threads(file_size, file_type, sort_type):
    """
    Mainly for adjusting the os_mem and os_threads for MacOS which has a limit of 256 open files per process.
    samtools sort creates temporary files equal to 1/2.3 of the memory per thread.  But if that
    memory size will create more than 256 temporary files, sort will fail. So adjust memory per thread
    up and # of threads down, as needed, to keep the temporary files under 240. Multiplier is based on BAM file
    size.  Adjust if CRAM file size given.  Name sort is ~50% larger total temp files than coordinate sort.

    Remember: file_size and os_totmem are in bytes; os_mem is in millions of bytes and a string
    """

    DEBUG(f'In Adjust_Mem for Sort with File Size: {file_size // 10 ** 9} GB, {file_type}, "{sort_type}" Sort')

    # Constant values / multipliers to help adjust values
    temp_mult = 2.3         # If memory per thread is x, then actual temp file size per thread is x/2.3
    cram_mult = 1.95 if file_type == "CRAM" else 1       # Adjust file size to likely BAM file size if CRAM to start
    name_mult = 1.45 if sort_type == "Name" else 1       # Name sorting needs more temp file space than coordinate
    files_max = 256 - 10 + wgse.os_threads               # For MacOS, max 256 open files per process. Less as some
    #                                   opened by process. More as # threads mem blocks kept out of temp file space

    # Calculate minimum mem per thread on MacOS to avoid exceeding 250 maximum open files during sort
    min_mpt = int(file_size * temp_mult * cram_mult * name_mult / files_max)  # Minimum mem per thread on MacOS
    if min_mpt == 0:      # Cannot adjust parameters if an (near) empty BAM file is being used to sort
        return "0M", 0

    # Let user know how much free file space we need; giving them a chance to make it available or cancel out
    # app = f'{sort_type} Sort'
    # isize = file_size * cram_mult / 10**9               # Nominal size of final BAM (in GB)
    # tempsize = str(int(isize * name_mult + isize))      # Size of raw BAM from aligner + coord sorted BAM pre-final (GB)
    # outsize = str(int(isize))                           # Size of final BAM (assume BAM always for worst case) (GB)
    # message = wgse.lang.i18n["infoFreeSpace"].replace(
    #    "{{APP}}", app).replace("{{SIZE}}", tempsize).replace("{{FINAL}}", outsize)
    # if not wgse_message("okcancel", "infoFreeSpaceTitle", True, message):
    #    return "0M", 0

    cur_mpt = wgse.os_totmem // wgse.os_threads     # note: os_mem is a string and in millions; this is in bytes

    if wgse.os_plat == "Darwin" and min_mpt > cur_mpt:

        # We need to adjust the memory per thread (mpt) up and # of threads down to prevent too many temp files in sort
        os_threads = wgse.os_totmem // min_mpt      # Adjust OS threads downward to allow increase of mem per thread

        # If 0 threads, cannot increase mem per thread to the needed value so print error and error out
        if os_threads == 0:
            size = str(round(file_size * cram_mult / 10 ** 9, 0))
            wgse_message("error", "errNoMemMacOSTitle", True, wgse.lang.i18n["errNoMemMacOS"].replace("{{SIZE}}", size))
            return "0M", 0

        # We can adjust Mem per thread further upward due to the small, discrete number of processors
        min_mpt = wgse.os_totmem // os_threads

        os_mem = str((min_mpt // 10 ** 8) * 100 if min_mpt > 10 ** 8 else (min_mpt // 10 ** 6)) + 'M'

        DEBUG(f'*** Adjust "{sort_type}" sort command on MacOS to use memory per thread: {wgse.os_mem} to {os_mem}, '
              f'(Perf) Threads: {wgse.os_threads} to {os_threads}')
        return os_mem, os_threads

    else:       # No adjustment necessary; use values set by platform at program startup or override by user earlier
        DEBUG(f'*** No adjustment. Using memory per thread: {wgse.os_mem} and (Perf) Threads: {wgse.os_threads}')
        return wgse.os_mem, wgse.os_threads


def button_sort_BAM():
    """
        Recreate BAM in coordinate sorted format. This is an exception where generated data is put where the current
        BAM file is; not in the output or temp directory.  We will use the sorted BAM going forward.
        Button is not available to the user unless the BAM is unsorted.
    """

    if wgse.BAM.Sorted is None:         # Unaligned BAM cannot have stats or sort commands run on it
        wgse_message("error", 'errUnalignedBAMTitle', False, 'errUnalignedBAM')
        return

    samtools = wgse.samtoolsx_qFN
    bamfile  = wgse.BAM.file_qFN
    tempdir  = f'"{wgse.tempf.FP}"'

    # Try to avoid adding conflicting names
    out_FPB =  wgse.BAM.file_FPB.replace("_unsorted", "_sorted") if "_unsorted" in wgse.BAM.file_FPB else \
            f'{wgse.BAM.file_FPB}_sorted'

    # CRAM decode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    suffix = ".cram" if wgse.BAM.file_type == "CRAM" else ".bam"
    out_qFN = f'"{out_FPB}{suffix}"'

    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist

    sort_mem, sort_cpus = _adjust_mem_threads(wgse.BAM.file_stats.st_size, wgse.BAM.file_type, "Coord")
    if sort_cpus == 0:         # Not enough memory to run samtools sort on MacOS; already reported the error
        mainwindow_resume()
        return

    # Samtools sort cannot accept a CRAM and have a reference genome specified so view first and pipe to sort
    # Although not clear if --reference works on sort command; does not complain but does not seem to help
    #   samtools view | sort on CRAM is 374min (67min real) on i9 - 9900 with 40GB DRAM and NVMe SSD
    #   samtools sort --reference on CRAM directly is 434m (119min real)
    #   sort on BAM is 370min (45min real) (BAM made from CRAM earlier)
    commands = (
        f'{samtools} view -uh --no-PG {cram_opt} {bamfile} | '
        f'  {samtools} sort -T {tempdir} -m {sort_mem} -@ {sort_cpus} -o {out_qFN} \n'
    )

    run_bash_script("GenSortedBAM", commands)

    set_BAM_file(unquote(out_qFN))     # Replace selected BAM and garbage collect old one

    mainwindow_resume()


def button_index_BAM():
    """
        Create BAM Index file. This is an exception where generated data is put where the current
        BAM file is; not in the output or temp directory.
    """
    # This is an exception where generated data is put where the BAM file is; not the Output or Temp directories
    samtools = wgse.samtoolsx_qFN
    bamfile  = wgse.BAM.file_qFN
    commands = f'{samtools} index {bamfile}\n'

    run_bash_script("GenBAMIndex", commands)

    if wgse.BAM.check_for_bam_index():  # if successful
        wgse.BAM.Indexed = True
        # For BAMs, run Stats (internal only) now that Index is available (saves clicking Stats button; only 1 sec)
        try:
            if wgse.BAM.Sorted is not None:     # Cannot run stats on an Unaligned BAM
                wgse.BAM.get_samtools_idxstats(button_directly=False)
        except:  # Error generating idxstats file when tried
            pass    # No need to report issue. Just trying an internal convenience run
        set_BAM_window_settings()  # If index became available, change its buttons state. As well as any stats filled in
    else:
        wgse_message("warning", 'InvalidBAMTitle', False, 'InfoBamIndexError')

    mainwindow_resume()


def button_CRAM_to_BAM():
    """
    Routine to convert a CRAM file to a BAM.  Both same SAM content; just different compression format.
    Internal call so only made if currently loaded file is a CRAM (vs a BAM).
    """
    # CRAM decode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist

    CRAM_qFN = wgse.BAM.file_qFN
    BAM_FN = f'{wgse.outdir.FPB}.bam'
    BAM_qFN = f'"{BAM_FN}"'
    DEBUG(f"Output BAM (final): {BAM_FN}")

    samtools = wgse.samtoolsx_qFN       # Know we are using 64 bit samtools 1.10; 64 bit required for CRAM decode
    commands = (
        f'{samtools} view -bh {cram_opt} -@ {wgse.os_threads} -o {BAM_qFN} {CRAM_qFN} \n'
        f'{samtools} index {BAM_qFN} \n'
    )

    run_bash_script("CRAMtoBAM", commands)

    set_BAM_file(BAM_FN)    # Replace selected BAM and garbage collect old one

    mainwindow_resume()


def button_BAM_to_CRAM():
    """
    Routine to convert a BAM file to a CRAM.  Both same SAM content; just different compression format.
    Internal call so only made if currently loaded file is a BAM (vs a CRAM).
    """
    # CRAM encode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}'
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist

    BAM_qFN = f'"{wgse.BAM.file_FN}"'
    CRAM_FN = f'{wgse.outdir.FPB}.cram'
    CRAM_qFN = f'"{CRAM_FN}"'
    DEBUG(f"Output CRAM (final): {CRAM_FN}")

    samtools = wgse.samtoolsx_qFN  # Know we are using 64 bit samtools 1.10; 64 bit required for CRAM encode
    commands = (
        f'{samtools} view -Ch {cram_opt} -@ {wgse.os_threads} -o {CRAM_qFN} {BAM_qFN} \n'
        f'{samtools} index {CRAM_qFN} \n'
    )

    run_bash_script("BAMtoCRAM", commands)

    set_BAM_file(CRAM_FN)   # Replace wuth new CRAM and garbage collect old one

    mainwindow_resume()


def button_realign_BAM(paired_BAM=True):
    """
        Button to ask to realign BAM/CRAM/SAM from one reference model to another. Involves doing an unalign
        back to FASTQ files and then alignment to a new BAM/CRAM/SAM after automatically selecting a "liftover"
        reference model.  This button is fixed to do to/from Build 37/38 only of the same class model.
        All functionality pushed into the button_align_command..
    """
    def _realign_exit(success=False):
        wgse.BAM = save_BAM if not success else ""     # Already posted message before call
        mainwindow_resume()
        return success

    save_BAM = wgse.BAM     # Will save current (if it exists) and let it garbage collect on exit if replaced

    # Run ask align parameters now simply to confirm what will be used. This way, is at the start and the
    # computer can be left alone to complete this job without interruption. Will update info in current BAM
    # but wipe out if fails when we restore the saved BAM.
    confirmed, fastqs, new_refgenome, newBAM_oFN = ask_align_params(inRealign=True, inAlign=False)
    if not confirmed:
        return _realign_exit(success=False)

    # Todo need to pass along possibly modified values (not just confirmed)?

    time_expected = str(5+160/wgse.os_threads)
    wgse_message("warning", "RealignBAMTimeWarnTitle", True,
                 wgse.lang.i18n["RealignBAMTimeWarnMesg"].replace("{{time}}", time_expected))

    # Need FASTQs (unaligned) so simply make button call. Returns FASTQs in wgse.BAM.Rxfastq variables
    # Will check if already exist and return FALSE for made_new if so.
    fastqs_exist, made_new_fastqs = button_unalign_BAM(inRealign=True)
    if not fastqs_exist:
        return _realign_exit(success=False)

    # This is an 8 hour to 6 day job; split into parts in the routine
    made_new_BAM = button_align_BAM(inRealign=True)  # Will automatically load newly created BAM/CRAM

    # If successfully replaced BAM and made new FASTQs as part of it; then delete FASTQs
    if made_new_fastqs and made_new_BAM and wgse.BAM and save_BAM and wgse.BAM != save_BAM:
        wgse.tempf.list.append(nativeOS(save_BAM.R1fastq_FN));  save_BAM.R1fastq_FN = None
        wgse.tempf.list.append(nativeOS(save_BAM.R2fastq_FN));  save_BAM.R2fastq_FN = None
        save_BAM = None     # Garbage collect old BAM class structure

    return _realign_exit(success=True)


def enough_free_space_align(fastq_size, newBAM_FS="BAM", fastq_exist=False, report=False):
    """ FASTQ is total size of both, if paired-end. """
    # If fastqs do not yet exist, then they are about to be created; so add their size to the outneed
    mult = 2 if fastq_exist else 3

    bam_size = int(fastq_size)                                      # Nominal size of final BAM
    cram_size = int(bam_size / 2 if newBAM_FS == ".cram" else 0)    # Additional size if creating a CRAM after
    tempneed = int(round( bam_size * 6 / 10**9, 0))    # Coordinate sort needs fully uncompressed BAM in temp file area
    outneed = int(round( (bam_size * mult + cram_size) / 10**9, 0))    # Size of RAW, final BAM, & final CRAM

    ofree = int(round(wgse.outdir.update_free() / 10**6, 0)) if wgse.outdir.FP else 10**7      # 10 TB if not set yet
    tfree = int(round(wgse.tempf.update_free() / 10**6, 0))
    samedisk = tfree == ofree       # Get rough MB to check if the same; will need many GB

    ofree = int(round(ofree / 10**3, 0))    # Convert to GB
    tfree = int(round(tfree / 10**3, 0))

    ostr = f'{outneed:,} GB'
    tstr = f'{tempneed:,} GB'

    DEBUG(f"Alignment Free Disk Check: Temp and Outdir same? {samedisk}; \n"
          f"\t\ttemp dir: {tfree:,} GB avail vs {tstr} needed; \n"
          f"\t\tout dir: {ofree:,} GB avail vs {ostr} needed")

    enough = (tempneed + outneed) < tfree if samedisk else tempneed < tfree and outneed < ofree
    if report and not enough:
        outneed = f'{outneed} GB'              # Converted to GB string
        tempneed = f'{tempneed} GB'            # Converted to GB string
        message = wgse.lang.i18n["infoFreeSpace"].replace(
            "{{APP}}", "Align").replace("{{SIZE}}", tstr).replace("{{FINAL}}", ostr)
        if wgse_message("okcancel", "infoFreeSpaceTitle", True, message):
            enough = True   # User told us to proceed anyway; or maybe made the space available (no recheck)

    return enough, tempneed, outneed


def ask_align_params(inRealign=False, inAlign=False):
    """
    GUI pop-up to set the parameters for an alignment. As well as indicate if there is enough free space.
    Will be called even from the Realign command (but pre-filled with all the values). When a FASTQ file
    name is selected, and if the Output Directory is not yet selected, a default outdir will be created and used
    just as happens when loading a BAM file. Any BAM File with a path is stripped of the path. The BAM file will
    always be in the output directory. If a fastq and refgenome are selected before a bamfile by the user,
    a default bamfile name (with .cram extension) is constructed and displayed to make it easier for the user.
    """

    # Local function to enable setting of local variables; should we change to a class?
    def setting_click(action, **args):
        nonlocal fastq1, fastq2, refgen, bamfile, bamset, fastq_good

        # Do not allow values to be changed if inRealign
        result = action(**args)     # Execute action specified with args

        if action.__name__ == "ask_fastq_to_align":
            fastq1 = result[0] if len(result) > 0 else ""
            fastq2 = result[1] if len(result) > 1 else ""
            fastq_good = update_fastq_labels()

            # Like for BAM, setup a default outdir if not yet set based on the fastq1 file location
            if wgse.outdir and not wgse.outdir.user_set and fastq1:
                wgse.outdir.new_default(os.path.dirname(fastq1))

            update_needed_space() if fastq1 else ""

            # if no bamfile set, then set a default bamfile for the user (they can always override; default is cram)
            if fastq1 and wgse.outdir and wgse.outdir.FP and not bamset:
                bambase = _name_from_FASTQ(os.path.basename(fastq1))
                bamfile = f'{wgse.outdir.FP}{bambase}.cram'
                if refgen:
                    bamfile = mod_bamfile_with_refgen(bamfile)
                newbamSelectLabel.configure(text=os.path.basename(bamfile))
                # Still default so do not set bamset

        elif action.__name__ == "ask_refgenome":
            # Note: do not require the refgen file be available; will query the user to download later
            refgen = result
            fgc = maroon if not wgse.reflib.get_refgenome_qFN(refgen) else wgse.deflab_fg
            refgenSelectedLabel.configure(text=refgen, fg=fgc)

            # If BAM not set by user and a default value; then modify the default for the new refgenome selected
            if not bamset and bamfile:
                bamfile = mod_bamfile_with_refgen(bamfile)
                newbamSelectLabel.configure(text=os.path.basename(bamfile))

        elif action.__name__ == "ask_BAM_filename":
            bamfile = result
            base = os.path.basename(bamfile)
            base, ext = os.path.splitext(base)
            newbamSelectLabel.configure(text=base)
            bamset = True       # Once set by a user; no longer update the "default" value

            # If no outdir set yet, try to set a default based on the fastq1 file location
            if fastq1 and wgse.outdir and not (wgse.outdir.user_set or wgse.outdir.FP):
                wgse.outdir.new_default(os.path.dirname(fastq1))

            # Even if user set a path for the bamfile, we want to override and put it in the outdir
            if wgse.outdir and wgse.outdir.FP:
                bamfile = f'{wgse.outdir.FP}{base}'    # Put new BAM file in an existing outdir

            if ext == ".CRAM":
                update_needed_space()       # Need increases if user specified a CRAM

        # Only enable the confirm button once everything is set correctly; key being fastq's exist for the align
        if fastq1 and refgen and bamfile and (fastq_good if not inRealign else True):
            confirmedButton.configure(state="normal")

        askAlignWindow.update()

    def ask_confirm():
        nonlocal confirmed
        # We only make the confirmed button available when everything is OK, so just close window
        try:
            askAlignWindow.destroy()
        except:
            pass
        confirmed = True

    def ask_cancel():
        try:
            askAlignWindow.destroy()
        except:
            pass

    def check_fastq(fastq):
        return fastq and os.path.exists(fastq) and os.path.isfile(fastq) and os.path.getsize(fastq) > 200 * 10*6

    def update_fastq_labels():
        nonlocal fastq1, fastq2

        fq1_good = False
        fq2_good = False

        if fastq1:
            fq1_good = check_fastq(fastq1)
            color = red if not inRealign else maroon      # if inRealign, fastq not required (so pink)
            fg1 = color if not fq1_good else wgse.deflab_fg
            fastqSelectedLabel1.configure(text=os.path.basename(fastq1), fg=fg1)

        if fastq2:
            fq2_good = check_fastq(fastq2)
            color = red if not inRealign else maroon      # if inRealign, fastq not required (so pink)
            fg2 = color if not fq2_good else wgse.deflab_fg
            fastqSelectedLabel2.configure(text=os.path.basename(fastq2), fg=fg2)

        return fastq1 and fq1_good and (fq2_good if fastq2 else False)

    def update_needed_space():
        bamfile_FS = os.path.splitext(os.path.basename(bamfile))[1] if bamfile else "BAM"

        # Set the fastq size as all additional needed space is based on that.
        fsize = (os.path.getsize(fastq1) if fastq1 and os.path.exists(fastq1) and os.path.isfile(fastq1) else 0) + \
                (os.path.getsize(fastq2) if fastq2 and os.path.exists(fastq2) and os.path.isfile(fastq2) else 0)

        # If inRealign and FASTQs have not yet been created, estimate fsize from the curent BAM size
        fastq_exist = fsize > 0
        fsize = wgse.BAM.file_stats.st_size if inRealign and wgse.BAM and fsize == 0 else fsize

        enough, tneed, oneed = enough_free_space_align(fsize, bamfile_FS, fastq_exist, report=False)
        needfg = maroon if not enough else wgse.deflab_fg
        tneedLabel.configure(text=f'Temp: {tneed} GB', fg=needfg)
        oneedLabel.configure(text=f'OutDir: {oneed} GB', fg=needfg)

        return enough

    def mod_bamfile_with_refgen(bfile):
        # Modify the default bamfile value to have a new refgen code (or replace one that may be there)
        found = False
        codes = wgse.reflib.refgenome_codes()
        for code in codes:
            if f'_{code}' in bfile:  # Do even if refgen = code so sets found=true
                bfile.replace(code, refgen)
                found = True
                break       # Only replaces one hyphenated _refgen code in the original fastq name; not if 2 or more

        # If not found, then lets add the refgen code to the name so we do not clobber an existing with the default
        if not found:
            base, ext = os.path.splitext(bfile)
            bfile = f'{base}_{refgen}{ext}'

        return bfile

    # Set initial / seed values based on loaded BAM (if any)
    if wgse.BAM:
        if not wgse.BAM.R1fastq_FN:
            wgse.BAM.find_and_set_FASTQs(inRealign=True)        # If in realign, always set BAM.Rxfastq_FN on return
        fastq1, fastq2 = [wgse.BAM.R1fastq_FN, wgse.BAM.R2fastq_FN]
    else:
        fastq1, fastq2 = ["", ""]
    refgen, _ = wgse.reflib.refgen_liftover(wgse.BAM.Refgenome) if wgse.BAM and wgse.BAM.Refgenome else ("", 0)
    bamfile = wgse.BAM.realign_BAM_filename(refgen) if refgen else ""

    # If ..., then return default parameters as they were confirmed in the Realign call to this routine
    if inRealign and inAlign:
        return True, [fastq1, fastq2], refgen, bamfile

    bamset = True if bamfile else False
    confirmed = False       # Set when confirmed button hit and all is good. Button only enabled once everything is set.
    fastq_good = False      # Purely for not inRealign as FASTQs must exist, be large enough, etc.
    free_good = False       # TO indicate if not enough free space; but use as advisory for now

    # ---------- Start of main window creation to Confirm and possibly adjust values ---------------------------------

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    askAlignWindow = Toplevel(wgse.window)
    askAlignWindow.transient(wgse.window)
    askAlignWindow.title(wgse.lang.i18n['AskAlignTitle'])
    askAlignWindow.geometry(offset(wgse.window.geometry(), 50))
    askAlignWindow.protocol("WM_DELETE_WINDOW", ask_cancel)
    askAlignWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    askAlignWindow.columnconfigure(0, weight=1)
    askAlignWindow.rowconfigure(0, weight=1)

    Label(askAlignWindow, text=f"{wgse.lang.i18n['AskAlign']}",
          font=font['14b']).grid(row=0, column=0, columnspan=3, padx=5, pady=2); crow = 1

    fastqSelectButton = Button(askAlignWindow, text=wgse.lang.i18n['SelectFASTQFiles'], font=font['14'],
                               command=partial(setting_click, ask_fastq_to_align))
    fastqSelectButton.grid(column=0, row=crow, padx=5, pady=2, rowspan=2, sticky=W)
    fastqSelectedLabel1 = Label(askAlignWindow, text="", font=font['14'])
    fastqSelectedLabel1.grid(column=1, row=crow, padx=5, pady=2, columnspan=2); crow += 1
    fastqSelectedLabel2 = Label(askAlignWindow, text="", font=font['14'])
    fastqSelectedLabel2.grid(column=1, row=crow, padx=5, pady=2, columnspan=2); crow += 1

    fastq_good = update_fastq_labels()

    fgc = maroon if not wgse.reflib.installed_refgenome(ref_genome=refgen) else wgse.deflab_fg
    refgenSelectButton = Button(askAlignWindow, text=wgse.lang.i18n['SelectReferenceGenome'], font=font['14'],
                                command=partial(setting_click, ask_refgenome, inBAM=False, window=askAlignWindow))
    refgenSelectButton.grid(column=0, row=crow, padx=5, pady=2,  sticky=W)
    refgenSelectedLabel = Label(askAlignWindow, text=refgen, font=font['14'], fg=fgc)
    refgenSelectedLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=2); crow += 1

    bam_base = os.path.basename(bamfile) if bamfile else ""

    newbamSelectButton = Button(askAlignWindow, text=wgse.lang.i18n['SelectBamFile'], font=font['14'],
                                command=partial(setting_click, ask_BAM_filename, new=True))
    newbamSelectButton.grid(column=0, row=crow, padx=5, pady=2, sticky=W)
    newbamSelectLabel = Label(askAlignWindow, text=bam_base, font=font['14'])
    newbamSelectLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=2); crow += 1

    Label(askAlignWindow, text=wgse.lang.i18n['DiskSpaceNeeded'], font=font['14']
          ).grid( column=0, row=crow, padx=5, pady=2)
    tneedLabel = Label(askAlignWindow, text=f'Temp: 0 GB', font=font['14'], fg=wgse.deflab_fg)
    tneedLabel.grid(column=1, row=crow, padx=5, pady=2,)
    oneedLabel = Label(askAlignWindow, text=f'OutDir: 0 GB', font=font['14'], fg=wgse.deflab_fg)
    oneedLabel.grid(column=2, row=crow, padx=5, pady=2); crow += 1

    free_good = update_needed_space()        # if inRealign, know the needed sizes already

    Button(askAlignWindow, text=wgse.lang.i18n['Cancel'], font=font['14'],
           command=ask_cancel).grid(row=crow, column=1, padx=5, pady=2)
    confirmedButton = Button(askAlignWindow, text=wgse.lang.i18n['Confirm'], font=font['14'],
                             state="disabled", command=ask_confirm)
    confirmedButton.grid(row=crow, column=2, padx=5, pady=2); crow += 1

    # Enable confirmed button if all three initial values have a default
    if bamfile and refgen and fastq1 and (fastq_good if not inRealign else True):
        confirmedButton.configure(state="normal")

    # Do not allow values to be changed in a realign; simply confirm
    if inRealign:
        fastqSelectButton.configure(state="disabled")
        refgenSelectButton.configure(state="disabled")
        newbamSelectButton.configure(state="disabled")

    # Display and wait for Confirm or Close button
    askAlignWindow.update()
    askAlignWindow.grab_set()          # Toplevel equivalent of mainloop
    askAlignWindow.wait_window()

    return confirmed, [fastq1, fastq2], refgen, bamfile


def button_align_BAM(inRealign=False):
    """
        Button to align (create) BAM/CRAM/SAM from a reference model and FASTQs.
        If button on GUI, then generate pop-up(s) asking for values (FASTQs, target file name and ref genome)
        Otherwise, if part of a Realign for the current BAM, then can determine all needed values here.
        Knows to use minimap2 on long reads (Oxford Nanopore, single Fastq). If single-read, then 2nd
        FASTQ file variable will be set null.

        Major Steps once parameters are setup (each with its own PleaseWait):
          (a) Index Reference Genome for alignment (if not already) (1+ hour, ~5GB)
          (b) Call ALignment tool on FASTQs to create initial raw BAM (8+ hours, ~50GB)
          (c) Cleanup (fixmate, markdups, sort) and Index new BAM (1 hour, same size, if short-read)
          (d) If needed, convert to CRAM and remove BAM (1 hour)

        Inspirations for Samtools bioinformatics pipeline used here:
        https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html
        https://www.htslib.org/workflow/ (first section on FASTQ to BAM/CRAM)
        https://gist.github.com/tkrahn/7dfc51c2bb97a6d654378a21ea0a96d4 (although some he copied from us :)
    """

    # Todo should main content be moved to bamfile.py?

    def _align_exit(success=False):
        wgse.BAM = saveBAM if not success else ""                               # Restore BAM in case removed
        wgse.outdir.change(saveOutFP, user_set=False) if not success else ""    # Restore (default) Outdir if changed
        mainwindow_resume() if wgse.window and wgse.dnaImage else ""
        return success

    # CAUTION: Only use wgse.BAM if in reAlign; otherwise may be set but not valid for direct Align button click
    saveBAM = wgse.BAM              # Save BAM to restore on error
    saveOutFP = wgse.outdir.FP      # Save outdir FP to restore on error; Outdir may change to FASTQs path if default

    # ---------------------------------------------------------------------------------------------------------------
    # Get and confirm key parameters from the user in a pop-up; if in Realign, do not report an issue
    confirmed, fastqs, new_refgenome, newBAM_oFN = ask_align_params(inRealign, inAlign=True)
    if not confirmed:
        return _align_exit()

    # ---- (1) Process FASTQ file and Related settings --------------------------------------------------
    # note: if one or both set then ...
    if fastqs:
        paired = len(fastqs) == 2 and fastqs[1]
        f1_FN = fastqs[0]
        f2_FN = fastqs[1] if paired else ""

    # else null out variables; will pickup and error out later (should not happen as user cannot clear defaults)
    else:
        paired = False
        f1_FN = f2_FN = ""

    DEBUG(f'FASTQ(s) to align: paired? {paired}, {os.path.basename(f1_FN)}, {os.path.basename(f2_FN)}')

    # Cheat and use BAM stats for needed FASTQ stats (if available)
    if inRealign and saveBAM and saveBAM.Stats:
        sequencer = saveBAM.Sequencer
        numsegs_flt = saveBAM.raw_segs_read         # Stored as float as may be greater than 4 GB
        read_length = saveBAM.avg_read_length

    # Process a FASTQ to get stats; only need one as if paired can double the number of segments
    else:
        sequencer, numsegs_flt, read_length = process_FASTQ(f1_FN, paired)

    numsegs = int(round(numsegs_flt, 0))
    DEBUG(f'Sequencer: {sequencer}, number of segments: {numsegs/10**6} M, read length: {read_length}')

    # Setup quoted FASTQ files; note special universal version for WSL needs if using WSL BWA
    f1_oFN = nativeOS(f1_FN)
    f1_quFN = f'"{universalOS(f1_FN, wgse.wsl_bwa_patch if wgse.os_plat == "Windows" else False)}"'
    f2_oFN = nativeOS(f2_FN) if paired else ""
    f2_quFN = f'"{universalOS(f2_FN, wgse.wsl_bwa_patch if wgse.os_plat == "Windows" else False)}"' if paired else ""

    # Finally check that FASTQs are large enough
    fastq_size = (os.path.getsize(f1_oFN) if os.path.isfile(f1_oFN) else 0) + \
                 (os.path.getsize(f2_oFN) if paired and os.path.isfile(f1_oFN) else 0)
    if fastq_size < 200 * 10**6:
        wgse_message("error", 'errFASTQMissingTitle', True,
                     f'{wgse.lang.i18n["errFASTQMissing"]}'.replace("{{f1}}", f1_FN).replace("{{f2}}", f2_FN))
        return _align_exit()

    # ---- (2) Process Reference Genome and related settings ----------------------------------------
    if not new_refgenome or new_refgenome in ["error", "unknown"]:   # Must have refgenome specified; so simply exit
        wgse_message("error", 'errRefGenNameMissingTitle', False, "errRefGenNameMissing")
        return _align_exit()

    # Check if reference genome file exists. If not, retrieve if they approve; or pause to give them a chance; or exit
    refgen_qFN = wgse.reflib.get_refgenome_qFN(new_refgenome)
    if wgse.reflib.missing_refgenome(refgen_qFN):
        # Above reports an error if the reference file does not exist
        return _align_exit()

    refgen_oFN = nativeOS(unquote(refgen_qFN))  # Remove quotes and change to native OS
    DEBUG(f'Reference Genome: {new_refgenome}, File: {os.path.basename(refgen_oFN)}')

    aligner = "minimap2" if "Nanopore" in sequencer else \
              "hisat2"   if "GRCh"     in new_refgenome else \
              "bwa"

    # --- (3) Set output BAM filename for new alignment -----------------------------------------------------
    if not newBAM_oFN:
        # Todo report error pop-up before returning here
        return _align_exit()

    # Define the final BAM / CRAM file name
    newBAM_FBS = os.path.basename(newBAM_oFN)
    newBAM_FB, newBAM_FS = os.path.splitext(newBAM_FBS)
    newBAM_FPB = f'{wgse.outdir.FP}{newBAM_FB}'
    newBAM_FN = f'{wgse.outdir.FP}{newBAM_FBS}'
    newBAM_qFN = f'"{newBAM_FN}"'
    DEBUG(f'New BAM File Basename: {newBAM_FB}, extension: {newBAM_FS}')

    # Define a trueBAM for alignment output in case newBAM is asking for a CRAM
    trueBAM_FN = newBAM_FN.replace(".cram", ".bam") if newBAM_FS == ".cram" else newBAM_FN
    trueBAM_qFN = f'"{trueBAM_FN}"'
    trueBAM_oFN = nativeOS(trueBAM_FN)

    # --- (4) Free File Space Requirements -------------------------------------------------------------------
    # Check if enough disk and memory free space to accomplish alignment, cleanup and sort
    cpus = wgse.os_threads
    # mem = wgse.os_mem

    enough, _, _ = enough_free_space_align(fastq_size, newBAM_FS, fastq_exist=True, report=True)
    if not enough:
        # Error reported in above if an issue; if returns false, then not enough space
        return _align_exit()

    # MacOS has a limitation in samtools sort due to only 255 open files per process. So use separate spec for sort.
    sort_mem, sort_cpus = _adjust_mem_threads(fastq_size, "BAM", "Coord")   # Reports error if an issue on MacOS
    if sort_cpus == 0:
        return _align_exit()

    # ---------------------------------------------------------------------------------------------------------------
    # Setup other file names (intermediate, reports, etc) and tool paths

    # Capture and save markdup report in output directory
    markdup_result_FN = f'"{newBAM_FPB}_markdup.txt"'

    # Important intermediate BAM files (moved to output directory to keep until very end; for crashes and hangs)
    outd_rawalign_FN  = f'{wgse.outdir.FP}{newBAM_FB}_raw.bam'
    outd_rawalign_qFN = f'"{outd_rawalign_FN}"'
    outd_rawalign_oFN = nativeOS(outd_rawalign_FN)
    # outd_rawalign_quFN = f'"{universalOS(outd_rawalign_FN,wsl_mode)}"'
    outd_sorted_qFN  = f'"{wgse.outdir.FP}{newBAM_FB}_sorted.bam"'
    outd_sorted_oFN = nativeOS(unquote(outd_sorted_qFN))

    # For temporary files area in sort and similar commands
    tempdir_qFN  = f'"{wgse.tempf.FP}"'

    # Special adjustments for WSL BWA mode
    wsl_mode = wgse.wsl_bwa_patch and wgse.os_plat == "Windows"
    if wsl_mode:     # if running WSL BWA, need WSL cmd paths
        bwa = 'wsl bwa'
        rg = '"@RG\\\\\\tID:1\\\\\\tSM:WGSE\\\\\\tLB:lb"'     # WSL takes a swipe at escape processing

        # bgzip = f'"{universalOS(unquote(wgse.bgzipx_qFN), wsl_mode)}"'
        # python = f'"python3"'
        # progress = f'"{universalOS(progress, True)}"'   # Change to WSL format
    else:
        bwa = wgse.bwax_qFN
        rg = '"@RG\\tID:1\\tSM:WGSE\\tLB:lb"'
        # Special file path handling in universalOS for WSL command line; if needed

    refgen_quFN = f'"{universalOS(unquote(refgen_qFN), wsl_mode)}"'

    # Simpified Tool references
    bgzip = wgse.bgzipx_qFN
    python = wgse.python3x_qFN
    progress = f'{wgse.prog_FP}progress.py'
    samtools = wgse.samtoolsx_qFN
    minimap2 = wgse.minimap2x_qFN

    prefix = f'/usr/bin/env TZ="" ' if wgse.os_plat == "Windows" else ""

    # --- (a) Alignment Reference Indices ------------------------------------------------------
    # Check for Alignment Index files -- reuse or generate if missing
    if aligner == "bwa" and os.path.isfile(refgen_oFN + ".bwt") and \
       os.path.getmtime(refgen_oFN + ".bwt") > os.path.getmtime(refgen_oFN):
        DEBUG(f'Using previous BWA Indices for RefGenome: {refgen_qFN}')

    elif aligner == "bwa":     # Generate BWA index
        # Index the Reference Genome (note: does not work if EBI Reference genome)
        # Todo Should we use "-a bwtsw" instead of default "-a is" ?
        debug_mode = "DEBUG" if wgse.DEBUG_MODE else ""
        commands = f'{bwa} index {refgen_quFN} 2> >({prefix} {python} {progress} - 0 {debug_mode} >&2 ) \n'
        run_bash_script("CreateAlignIndices", commands)

    elif aligner == "hisat2":
        # Todo If EBI reference model, need to create Hisat2 index?
        pass

    elif aligner == "minimap2":
        # Todo any prep for minimap2 aligner?  when to use / select
        pass

    # --- (b) Actual Alignment -----------------------------------------------------------------
    # 4 hours to 6+ days depending on the number of CPU cores
    # todo add picard and GATK markdup function in place of samtools (or as option)
    if os.path.isfile(outd_rawalign_oFN) and \
       os.path.getmtime(outd_rawalign_oFN) > os.path.getmtime(f1_oFN) and \
       0.8 < os.path.getsize(outd_rawalign_oFN) / fastq_size < 1.2:
        DEBUG(f'Found previous RAW alignment to reuse: {outd_rawalign_qFN}')  # Not likely but worth a try
    else:
        # Run the alignment command; compress the output simply because it is so large otherwise
        if aligner == "bwa":
            # If single read, then f2_quFN will be "" ; BWA understands if only a single FASTQ file paramter
            # Note: could be wsl bwa or native bwa; appropriate values should be setup ahead of time
            debug_mode = "DEBUG" if wgse.DEBUG_MODE else ""
            commands = (
                f'{bwa} mem -t {cpus} -R {rg} {refgen_quFN} {f1_quFN} {f2_quFN} '
                f'  2> >({prefix} {python} {progress} - {numsegs} {debug_mode} >&2 ) | '
                f'   {bgzip} -@ {cpus} > {outd_rawalign_qFN}\n'
            )
            #

        elif aligner == "minimap2":
            # For long read, go directly to making TrueBAM file (no cleanup step (c)). No ref index file either.
            if sort_cpus == 0:  # Not enough memory for sort; simply return as already reported error
                return _align_exit()

            commands = (
                f'{minimap2} -ax map-ont -t {cpus} {refgen_quFN} {f1_quFN} |'
                f'  {samtools} sort -T {tempdir_qFN} -m {sort_mem} -@ {sort_cpus} -o {trueBAM_qFN} \n'
                f'{samtools} index {trueBAM_qFN} \n'
            )

        else:   # if aligner == "hisat2":
            # Todo handle hisat2 alignment for GRCh / EBI models; really should error out
            commands = ""

        # Todo add pbmm aligner for PacBio HiFi CCS long-read FASTQ files; for reprocessing the T2T / HPP files
        run_bash_script("ButtonAlignBAM", commands)

    # --- (c1) Cleanup (fixmate, coordinate sort) ----------------------------------------------
    if os.path.isfile(outd_sorted_oFN) and os.path.isfile(outd_rawalign_oFN) and \
       os.path.getmtime(outd_sorted_oFN) > os.path.getmtime(outd_rawalign_oFN) and \
       0.5 < os.path.getsize(outd_rawalign_oFN) / os.path.getsize(outd_sorted_oFN) < 2:
        DEBUG(f'Found previous fixmate and sorted BAM to reuse: {outd_sorted_oFN}')

    elif aligner == "bwa":        # Cleanup built into minimap alignment stage
        # f'{samtools} view -uh --no-PG {outd_rawalign_qFN} | '     # Why did we have the view before Fixmate?
        if sort_cpus == 0:      # Not enough memory for sort; report error again as alignment puts out lots of messages
            fs_str = str(round(fastq_size / 10 ** 9, 0))
            wgse_message("error", "errNoMemMacOSTitle", True, wgse.lang.i18n["errNoMemMacOS"].replace("SIZE", fs_str))
            return _align_exit()

        commands = (
          f'{samtools} fixmate -m -O bam -@ {cpus} {outd_rawalign_qFN} - |'
          f'  {samtools} sort -T {tempdir_qFN} -m {sort_mem} -@ {sort_cpus} -o {outd_sorted_qFN} - \n'
        )
        run_bash_script("AlignCleanup", commands)
        # Note that a coordinate sorted BAM can be smaller than a name or unsorted BAM. Maybe
        #  compression is more effective when the similar sequences are brought closer to each other?

    # --- (c2) Cleanup2 (markdup and index new BAM) --------------------------------------------
    if os.path.isfile(trueBAM_oFN) and os.path.isfile(outd_sorted_oFN) and \
       os.path.getmtime(trueBAM_oFN) > os.path.getmtime(outd_sorted_oFN) and \
       0.5 < (os.path.getsize(trueBAM_oFN) / os.path.getsize(outd_sorted_oFN)) < 1.2:
        DEBUG(f'Found previous sorted, final BAM to reuse: {trueBAM_FN}')

    elif aligner == "bwa":        # Cleanup built into Long Read FASTQ files turned into BAMs
        # Check if samtools has read-coord option and sequencer has names that can be parsed for optical duplicates
        valid_markdup = any(elem in sequencer for elem in wgse.valid_markdup)
        if sequencer and valid_markdup:  # and (wgse.samtools_version[1] < 15 and "BGI" in sequencer):
            length = 100 if "HiSeq" in sequencer else 2500 if "Novaseq" in sequencer else 300  # for BGI DNB
            if wgse.samtools_version[1] > 14:
                read = wgse.sequencers[sequencer][1]
                order = wgse.sequencers[sequencer][2]
                coord_opts = f'--read-coords \"{read}\" --coords-order {order}'

            else:   # Old versions of samtools processed Illumina names with Row/Col specs embedded
                coord_opts = ""
                if "MGI" in sequencer:  # MGI files hang if try optical duplicate before version 15; simply turn off
                    length = 0

            markdup_opts = f'-d {length} {coord_opts}'
            # Todo set -l 400 (default is 300bp) to handle ySeq WG400?  What about nanopore? PacBio HiFi   ?

            commands = (
                f'{samtools} markdup -f {markdup_result_FN} {markdup_opts} '
                f'  -@ {cpus} -T {tempdir_qFN} {outd_sorted_qFN} {trueBAM_qFN} \n'
                f'{samtools} index {trueBAM_qFN} \n'
            )

        else:
            os.rename(outd_sorted_oFN, trueBAM_oFN)
            commands = f'{samtools} index {trueBAM_qFN} \n'

        run_bash_script("AlignCleanup2", commands)
        # Note that a coordinate sorted BAM is smaller than a name or unsorted BAM. Compression is more
        # effective when the similar sequences and quality strings are brought closer to each other

    # --- (d) Convert to CRAM and remove BAM (if a CRAM is requested) -------------------------
    if newBAM_FS == ".cram":
        # Note: would not find existing CRAM that is newer as we remove old BAM when converted
        # todo modify CRAM_to_BAM() to accept parameters and simply call that here
        commands = (
            f'{samtools} view -Ch -T {refgen_qFN} -@ {cpus} -o {newBAM_qFN} {trueBAM_qFN} \n'
            f'{samtools} index {newBAM_qFN} \n'
        )
        run_bash_script("BAMtoCRAM", commands)

        wgse.tempf.list.append(nativeOS(trueBAM_FN))            # Remove previous final BAM ...
        wgse.tempf.list.append(nativeOS(trueBAM_FN + ".bai"))   # and its index

    # else: No need to rename file as newBAM_FN and trueBAM_FN are identical names when final is a BAM

    # Replace selected BAM/CRAM and garbage collect old one; return success of change
    if not set_BAM_file(newBAM_FN):         # Will set wgse.BAM to newBAM and load; return sucess
        # Todo ? Delete raw and sorted BAMs if not in DEBUG?
        return _align_exit()

    # If finished successfully, then delete RAW and sorted bam file in output directory
    wgse.tempf.list.append(outd_rawalign_oFN)
    wgse.tempf.list.append(outd_sorted_oFN)
    wgse.BAM.R1fastq_FN = f1_FN     # Note: if in realign, and delete FASTQs, need to clear these
    wgse.BAM.R2fastq_FN = f2_FN

    # Todo check if Y-only BAM (effectively no / few reads on other primary sequences). If so, subset the BAM
    #  to just Y as reads may have aligned off the Y. Needed if started with Y only FASTQs or BAM

    return _align_exit(True)


def _name_from_FASTQ(fastq_file):
    """
    Try to be smart about naming reports on FASTQs.  Remove more than just the file extension if matching templates
    for common name extensions from delivered BAM files. Otherwise, simply strip common extension.
    """
    # Todo have similar code in the BAMFile class (also checks for existence); merge somehow?
    templates = ("_SA_L001_R1_001.fastq.gz", "_SA_L001_R2_001.fastq.gz",    # Dante Labs
                 "_1.fq.gz", "_2.fq.gz",                                    # Nebula Genomics
                 ".1.fq.gz", ".2.fq.gz",                                    # Sequencing
                 "_R1.fastq.gz", "_R2.fastq.gz",                            # WGSE unmap
                 ".fq.gz", ".fastq.gz")          # Simply remove the known extension, if any; maybe single-end?
    # (note: ySeq does not supply FASTQs)

    # matches = [x for x in templates if fastq_file.endswith(x)]        # list comprehension clearer?
    # return fastq_file if len(matches) == 0 else fastq_file[: -len(matches[0])]

    for template in templates:
        if fastq_file.endswith(template):
            return fastq_file[:-len(template)]      # Remove extension; leaving a BAM base name only

    # If get here; do not recognize FASTQ so simply return the original
    return fastq_file


def button_fastp_FASTQ():
    """
    Run the fastp command on existing FASTQ files; either attached to BAM or specified by the user.  Open the result.
    """

    # If BAM set and we can simply try and find it's matching FASTQs
    found = wgse.BAM and wgse.BAM.find_and_set_FASTQs()
    if found:     # Use found BAM FASTQ file(s)
        f1_FN = wgse.BAM.R1fastq_FN
        f2_FN = wgse.BAM.R2fastq_FN
        paired = wgse.BAM.ReadType == "Paired"
    else:                                       # or ask user for FASTQ(s) to use
        fastqs = ask_fastq_to_align()
        f1_FN = fastqs[0] if len(fastqs) >= 1 else ""
        f2_FN = fastqs[1] if len(fastqs) == 2 else ""
        paired = True if f2_FN else False

    f1_oFN = nativeOS(f1_FN)
    f2_oFN = nativeOS(f2_FN)

    if not (f1_FN and os.path.isfile(f1_oFN) and (not paired or f2_FN and os.path.isfile(f2_oFN))):
        wgse_message("error", 'errFASTQMissingTitle', True,
                     f'{wgse.lang.i18n["errFASTQMissing"]}'.replace("{{f1}}", f1_FN).replace("{{f2}}", f2_FN))
        mainwindow_resume()
        return

    fastp = wgse.fastpx_qFN
    if wgse.os_plat == "Darwin" and not fastp:
        wgse_message("error", 'errAppMissingTitle', True,
                     f'{wgse.lang.i18n["errAppMissing"]}'.replace("{{APP}}", "FASTP"))
        mainwindow_resume()
        return

    fastopt = f'-i "{f1_FN}"' + (f' -I "{f2_FN}"' if paired else "")

    # Name output file after BAM, if found, else strip FASTQ file name (smartly; not just split)
    f1_FB = wgse.BAM.file_FB if found else _name_from_FASTQ(os.path.basename(f1_FN))
    html = f'"{wgse.outdir.FP}{f1_FB}_fastp.html"'
    json = f'"{wgse.outdir.FP}{f1_FB}_fastp.json"'

    ohtml = nativeOS(unquote(html))
    if not os.path.isfile(ohtml):       # Only create file if it does not exist; takes 30-60 minutes so ...
        commands = f'{fastp} {fastopt} -h {html} -j {json} -R "{f1_FB} FASTP Report" \n'
        run_bash_script('ButtonFastp', commands)

    if os.path.isfile(ohtml):           # If file exists (whether just created or not), display it
        webbrowser.open_new(unquote(html))
    else:
        wgse_message("error", 'errFASTPMissingTitle', True,
                     f'{wgse.lang.i18n["errFASTPMissing"]}'.replace("{{FILE}}", html))
    mainwindow_resume()


def button_fastqc_FASTQ():
    """
    Run the fastqc java command on existing FASTQ files; either attached to BAM or specified by the user.
    Open the result in a webbrowser. If cannot find FASTQ files, ask for them. If output already created, simply
    open browser on result.
    """

    # If BAM set and we can find it's matching FASTQs, lets assume those are what they want.  Otherwise, ask which ones.
    found = wgse.BAM and wgse.BAM.find_and_set_FASTQs()
    if found:       # Use found BAM FASTQ file(s)
        f1_FN = wgse.BAM.R1fastq_FN
        f2_FN = wgse.BAM.R2fastq_FN
        paired = wgse.BAM.ReadType == "Paired"
    else:           # Else, ask user for FASTQ(s) to use
        fastqs = ask_fastq_to_align()
        f1_FN = fastqs[0] if len(fastqs) >= 1 else ""
        f2_FN = fastqs[1] if len(fastqs) == 2 else ""
        paired = True if f2_FN else False

    f1_oFN = nativeOS(f1_FN)
    f2_oFN = nativeOS(f2_FN)
    if not (f1_FN and os.path.isfile(f1_oFN) and (not paired or os.path.isfile(f2_oFN))):
        wgse_message("error", 'errFASTQMissingTitle', True,
                     f'{wgse.lang.i18n["errFASTQMissing"]}'.replace("{{f1}}", f1_FN).replace("{{f2}}", f2_FN))
        mainwindow_resume()
        return

    # Unlike fastp, paired files are processed independently. So no BAM naming of a single output file to start.
    ohtml1 = ohtml2 = None
    f1_FP, f1_FBS = os.path.split(f1_FN)
    f1_FB = fastq_base(f1_FBS)   # os.path.splitext(f1_FBS.replace(".fastq.gz", "").replace(".fq.gz", ""))[0]
    phtml1 = f'"{f1_FP}/{f1_FB}_fastqc.html"'
    pzip1  = f'"{f1_FP}/{f1_FB}_fastqc.zip"'
    fhtml1 = f'"{wgse.outdir.FP}{f1_FB}_fastqc.html"'
    fzip1  = f'"{wgse.outdir.FP}{f1_FB}_fastqc.zip"'
    ohtml1 = nativeOS(unquote(fhtml1))
    if paired:
        f2_FP, f2_FBS = os.path.split(f2_FN)
        f2_FB = fastq_base(f2_FBS)   # os.path.splitext(f2_FBS.replace(".fastq.gz", "").replace(".fq.gz", ""))[0]
        phtml2 = f'"{f2_FP}/{f2_FB}_fastqc.html"'
        pzip2  = f'"{f2_FP}/{f2_FB}_fastqc.zip"'
        fhtml2 = f'"{wgse.outdir.FP}{f2_FB}_fastqc.html"'
        fzip2  = f'"{wgse.outdir.FP}{f2_FB}_fastqc.zip"'
        ohtml2 = nativeOS(unquote(fhtml2))

    # Only run fastqc to create report if the report does not already exist; takes 30+ minutes to generate
    if not(os.path.isfile(ohtml1) and (not paired or os.path.isfile(ohtml2))):
        fastqc = wgse.fastqcx_qFN
        fastq_FN = f'"{f1_FN}"' + (f' "{f2_FN}"' if paired else "")

        if wgse.java17x_FN is None:  # Really still needed? Is not checked in Settings at startup?
            wgse_message("error", 'MissingJavaTitle', False, 'MissingJava')
            mainwindow_resume()
            return

        # Cannot get -D option in FastQC to work ...
        # fastopt = f'-Dfastqc.output_dir={wgse.outdir.oFP} -Djava.io.tmpdir={wgse.tempf.oFP} -Dfastqc.threads=2'
        commands = f'{wgse.java17x_FNp} {fastqc} {fastq_FN}\n'

        #  so have to move created files from FASTQ files location to Output Directory
        if f1_FP != wgse.outdir.FP:
            commands += f'{wgse.mvx_qFN} {phtml1} {fhtml1} ; {wgse.mvx_qFN} {pzip1} {fzip1}\n'

        # noinspection PyUnboundLocalVariable
        if paired and f2_FP != wgse.outdir.FP:
            # noinspection PyUnboundLocalVariable
            commands += f'{wgse.mvx_qFN} {phtml2} {fhtml2} ; {wgse.mvx_qFN} {pzip2} {fzip2}\n'

        run_bash_script('ButtonFastqc', commands)

    # If report still not available, then report error and return
    if not(os.path.isfile(ohtml1) and (not paired or os.path.isfile(ohtml2))):
        wgse_message("error", 'errFASTQCMissingTitle', True,
                     f'{wgse.lang.i18n["errFASTQCMissing"]}'.replace("{{FILE}}", fhtml1))
        mainwindow_resume()
        return

    # Special case for paired; use MultiQC (pip library routine) to merge report graphs into single output file
    if paired:
        # run MultiQC on two html / zip files to create single one; recreate even if already exists
        multiqc.run(wgse.outdir.FP, module=['fastqc'], outdir=nativeOS(wgse.outdir.FP), no_data_dir=True)

        # Like in Fastp above, be smart about naming combined report. Try and determine base, BAM file name from FASTQs
        f1_FB = _name_from_FASTQ(f1_FBS)
        html = f'{wgse.outdir.FP}{f1_FB}_multiqc.html'
        os.replace(f'{wgse.outdir.FP}multiqc_report.html', html)
    else:
        html = unquote(fhtml1)

    # If file exists (whether just created or not), display it
    if os.path.isfile(nativeOS(html)):
        webbrowser.open_new(html)
    mainwindow_resume()


def button_varqc_VCF():
    """
    Run the VariantQC java command (from DISCVRSeq) on the available FASTQ file(s).
    Note that DISCVRSeq requires jre8 (as does Picard and GATK3; which are embedded in it). Haplogrep and GATK4
    require jre11 or later.
    """

    # If BAM set and we can find it's matching VCFs, lets assume those are what they want.  Otherwise, ask which ones.
    found = wgse.BAM and wgse.BAM.find_VCFs()
    vcfs = wgse.BAM.VCFs if found else ask_vcfs_to_process()

    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "VarQC"))

    command = f'{wgse.java8x_FNp} "{wgse.jartools_FP}DISCVRSeq.jar" VariantQC '
    # Need the refgenome for (all) VCFs; if wgse.BAM, verify it is the same?
    # Need to create command line entry for each file; creating unique(?) "V" label
    # Need to create output file name; easy if wgse.BAM, ask user if not?

    # Needs oFP for all files passed to Java (C:\\... for Windows)
    # Needs .tbi index for any VCF: bcftools index -t *vcf.gz
    ''' Dante VCFs analysis:
    {wgse.java8x_FNp} "{wgse.jartools_oFP}DISCVRSeq.jar" VariantQC -R hs37d5*gz \
      -V:snp bamname.snp.vcf.gz -V:indel bamname.indel.vcf.gz -V:cnv bamname.cnv.vcf.gz \
      -V:sv bamname.sv.vcf.gz -O bamname_VarQC.html
    /mnt/c/wgse/dev/jre8/bin/java.exe -jar "C:\wgse\dev\jartools\DISCVRSeq.jar" VariantQC \
    -R "C:\wgse\dev\reference\genomes\hs37d5.fa.gz" \
    -V:snp1 608*Avanti.vcf.gz -V:snp2 608*snp.vcf.gz -O 60820188479441_VarQC.html
    '''
    ''' Nebula VCFs analysis:
    {wgse.java8x_FNp} "{wgse.jartools_oFP}DISCVRSeq.jar" VariantQC -R hs38*gz \
      -V:snp bamname.vcf.gz -O bamname_VarQC.html
    '''


def button_igv_BAM():
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "IGV"))
    # {javaw} --module-path={igv}/lib -Xmx8g -@{igv}/igv.args -Dproduction=true \
    #   -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true - \
    #   -module=org.igv/org.broad.igv.ui.Main --batch=batch.xml
    # batch.xml:
    #   new
    #   genome hg38
    #   load 10687_bwa-mem_hg38_sorted.bam
    #   goto chr1:21,100,100-21,100,150
    # Need to figure out how to specify our local indexed genome fasta files with annotations added;
    #   needs a json packaged file in gff it seems.
    # Biggest issue has been creating the "SNP - Gene name" / "search" track files for a WGS that is a reasonable size.


def button_unalign_BAM(inRealign=False):
    """
    Function to go from BAM / CRAM to FASTQ(s). Assumes BAM already subsetted if only want a subset FASTQ.
    Looks for result in output area first before acting
    """

    # Let's check if FASTQ's already exist; BAMfile class can do it and store locally if found
    missing_FASTQs = not wgse.BAM.find_and_set_FASTQs(inRealign)
    if missing_FASTQs:         # Could not find FASTQs (note: find_and_set leaves names in BAM class instance)

        r1fastq = r2fastq = sefastq = "/dev/null"
        if wgse.BAM.ReadType == "Paired":
            r1fastq = wgse.BAM.R1fastq_FN    # File names to create here
            r2fastq = wgse.BAM.R2fastq_FN
        else:
            sefastq = wgse.BAM.R1fastq_FN    # For single-end BAMs

        # CRAM file requires reference genome be specified with VIEW command
        cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
        samtools = wgse.samtoolsx_qFN
        bamfile  = wgse.BAM.file_qFN
        tempdir  = f'"{wgse.tempf.FP}"'

        if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
            mainwindow_resume() if not inRealign else ""
            return False, False      # Check routine reports error if reference does not exist

        sort_mem, sort_cpus = _adjust_mem_threads(wgse.BAM.file_stats.st_size, wgse.BAM.file_type, "Name")
        if sort_cpus == 0:  # Not enough memory to run samtools sort on MacOS; already reported the error
            mainwindow_resume() if not inRealign else ""
            return False, False

        # Sort in name order, then call fastq command to split and write FastQ's
        # Samtools sort cannot take the reference genome specification so have to view a CRAM first
        commands = (
            f'{samtools} view -uh --no-PG {cram_opt} {bamfile} | '
            f'  {samtools} sort -n -T {tempdir} -m {sort_mem} -@ {sort_cpus} -O sam | '
            f'  {samtools} fastq -1 "{r1fastq}" -2 "{r2fastq}" -0 "{sefastq}" -s /dev/null -n -@ {wgse.os_threads} \n'
        )

        run_bash_script('ButtonUnalignBAM', commands)

        # Todo error check on return

        if wgse.BAM.ReadType == "Paired":
            wgse.BAM.R1fastq_FN = unquote(r1fastq)
            wgse.BAM.R2fastq_FN = unquote(r2fastq)
        else:
            wgse.BAM.R1fastq_FN = unquote(sefastq)
            wgse.BAM.R2fastq_FN = ""

        mainwindow_resume() if not inRealign else ""
        return True, True       # FASTQs (now) available, FASTQs made as were missing

    mainwindow_resume() if not inRealign else ""
    return True, False          # FASTQs available, FASTQs reused as found


def button_select_outdir():
    """ Processing user button to specify Output Directory path """

    # Get directory to start in from the one containing the current set directory; if set. Otherwise, let OS decide
    initialdir = os.path.dirname(wgse.outdir.oFP[:-1]) if wgse.outdir and wgse.outdir.oFP else ""
    new_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectOutputDirectory'],
                                     initialdir=initialdir)      # Returns Unix/Universal, and no trailing slash
    if new_FP:      # Process user selection
        set_outdir(new_FP, user_set=True)
    else:
        DEBUG(f"Output Path not set (cancelled)")
    mainwindow_resume()


def button_unselect_outdir():
    if wgse.outdir:
        wgse.outdir.clear()                 # Just clear values; do not delete class
        outdirButton.configure(text=wgse.lang.i18n['Default'])
        outdirButton.configure(font=font['14'])
        outdirButton.configure(bg=wgse.defbut_bg)
        outdirUnselectButton.grid_remove()
        if wgse.BAM:                        # Have to reload to set default output directory (or should we just clear?)
            wgse_message("warning", "UnselectOutDirBAMTitle", False, "UnselectOutDirBAM")
            save = wgse.BAM  ;  wgse.BAM = None
            set_BAM_file(save.file_FN)      # If fails, just continue
        set_BAM_window_settings()           # Sets window labels of BAM file settings
        wgse.save_settings()                # Possibly changed BAM file and outdir; so save settings
    mainwindow_resume()


def set_outdir(new_FP, user_set=False):
    """"
        Set a new output directory (path) for all (permanent) output files.
        Called here from button_select_outdir, from standalone microarray module and the set_BAM_file
    """
    if not wgse.outdir or not new_FP or len(new_FP) < 3:
        DEBUG("*** Internal Error: Tried to set outdir before class setup, with an empty dir spec, or an invalid dir.")
        return

    if user_set:
        wgse.outdir.change(new_FP, user_set)      # Call class to actually change value
    else:
        wgse.outdir.new_default(new_FP)

    if wgse.outdir.FB and wgse.window and wgse.dnaImage:
        outdirButton.configure(text=wgse.outdir.FB if user_set else wgse.lang.i18n['Default'])
        outdirButton.configure(font=font['14b'] if user_set else font['14'])
        outdirButton.configure(bg=cyan if user_set else wgse.defbut_bg)
        outdirUnselectButton.grid()
    #else:
    #    wgse_message("error", 'InvalidPathWindowTitle', False, 'errOutputPathSpecialChars')

    wgse.save_settings()  # One of the saved settings; so go ahead and write out now (even if default so it is cleared)


def button_mtdna_BAM():
    """
    Processing user button to generate mito only BAM (always BAM, not CRAM)
    Put in output directory; does not replace currently selected BAM (extract button; not BAM settings button)
    """
    if wgse.BAM.Monly:    # MTonly BAM
        return True

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else  ""
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return False     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN
    outMTbam = f'"{wgse.outdir.FPB}_chrM.bam"'

    chromos = wgse.BAM.get_chr_name('M')

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {outMTbam} \n'

    run_bash_script("ButtonMitoBAM", commands)

    if not os.path.isfile(nativeOS(unquote(outMTbam))):
        wgse_message("error", 'errYTitle', False, 'errYBAMFile')

    mainwindow_resume()


def button_mtdna_VCF(internal=False):
    """
    Processing user button to generate the MITO VCF file (also from internal call to do the same).
    VCF is used by both mitoFASTA button and Haplogroup (haplogrep) call.  So reuse if available.
    Always place and leave in output directory for possible reuse.
    """
    if wgse.BAM.RefMito == "Yoruba" and not internal:
        wgse_message("warning", 'YorubaWindowTitle', False, 'YorubaWarning')

    chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
    chrM_oFN = nativeOS(unquote(chrM_qFN))

    if not (os.path.isfile(chrM_oFN) and os.path.getmtime(chrM_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        refgenome_qFN = wgse.BAM.Refgenome_qFN
        if wgse.reflib.missing_refgenome(refgenome_qFN):
            if not internal:
                mainwindow_resume()
            return False    # Check routine reports error if reference does not exist

        bcftools = f'{wgse.bcftoolsx_qFN}'
        tabix = f'{wgse.tabixx_qFN}'
        cpus = wgse.os_threads
        chrM = wgse.BAM.get_chr_name('M')
        ploidy = "1"  # mito-only

        #   Note -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
        commands = (
            f'{bcftools} mpileup -B -I -C 50 -r {chrM} -f {refgenome_qFN} -Ou {wgse.BAM.file_qFN} | '
            f'  {bcftools} call --ploidy {ploidy} -mv -P 0 --threads {cpus} -Oz -o {chrM_qFN} \n'
            f'{tabix} {chrM_qFN}\n'
        )

        run_bash_script("ButtonMitoVCF", commands)

        if not os.path.isfile(chrM_oFN):
            wgse_message("error", 'errMtTitle', False, 'errMtVCFFile')
            mainwindow_resume() if not internal else ""
            return False

    mainwindow_resume() if not internal else ""
    return True


def button_mtdna_FASTA():
    """
    Processing user button to generate the MITO FASTA file (also from internal call to do the same).
    Note: The FASTA in the industry is more like a FA and not a FASTQ. Not raw reads from the sequencer.
    But single ~16k base-pair segment of all calls. So like for Autosomal, have to mpileup and call.
    """
    if wgse.BAM.RefMito == "Yoruba":
        wgse_message("warning", 'YorubaWindowTitle', False, 'YorubaWarning')

    fasta_qFN = f'"{wgse.outdir.FPB}_mtdna.fasta"'
    fasta_oFN = nativeOS(unquote(fasta_qFN))
    if not (os.path.isfile(fasta_oFN) and os.path.getmtime(fasta_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
        if not button_mtdna_VCF(internal=True):
            mainwindow_resume()
            return      # Generate button reports error message if does not exist

        refgenome_qFN = wgse.BAM.Refgenome_qFN
        if wgse.reflib.missing_refgenome(refgenome_qFN):
            mainwindow_resume()
            return      # Check routine reports error if reference does not exist

        samtools = f'{wgse.samtoolsx_qFN}'
        bcftools = f'{wgse.bcftoolsx_qFN}'
        chrM = wgse.BAM.get_chr_name('M')

        # -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
        commands = f'{samtools} faidx {refgenome_qFN} {chrM} | {bcftools} consensus {chrM_qFN} -o {fasta_qFN} \n'

        run_bash_script("ButtonMitoFASTA", commands)

        if not os.path.isfile(fasta_oFN):
            wgse_message("error", 'errMtTitle', False, 'errMtFastaFile')

    mainwindow_resume()


def button_mtdna_haplogroup():
    """ Process user button to call haplogrep to calculate mtDNA haplogroup. """
    if wgse.BAM.RefMito == "Yoruba":
        wgse_message("error", 'YorubaWindowTitle', False, 'YorubaWarning')
        mainwindow_resume()
        return

    if wgse.java17x_FN is None:         # Really still needed? Is not checked in Settings at startup?
        wgse_message("error", 'MissingJavaTitle', False, 'MissingJava')
        mainwindow_resume()
        return

    chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
    if not button_mtdna_VCF(internal=True):
        mainwindow_resume()
        return  # button_mtdna_VCF will report an error message if it does not exist or is not created

    haplogrep = f'{wgse.java17x_FNp} "{wgse.jartools_FP}haplogrep.jar"'
    haplogroup_qFN = f'"{wgse.tempf.FP}mtdna_haplogroup.txt"'

    commands = f'{haplogrep} classify --in {chrM_qFN} --format vcf --out {haplogroup_qFN}\n'

    run_bash_script('ButtonMTHaplo', commands)

    haplogroup_oFN = nativeOS(unquote(haplogroup_qFN))
    if not os.path.isfile(haplogroup_oFN):
        wgse_message("error", 'errMtTitle', False, 'errMtHaploFile')
        mainwindow_resume()
        return

    mtdna_haplogroup = ""
    with open(haplogroup_oFN, "r") as source_content:
        for source_line in source_content:
            line_tabs = source_line.split("\t")
            if str(line_tabs[2]).strip('"') == "1":
                mtdna_haplogroup = line_tabs[1].strip('"')
                # Todo can we do a break here or might there be multiple with the third column equal to 1?
    report = wgse.lang.i18n['MitoHaploReport']\
        .replace("{{mthaplogroup}}", mtdna_haplogroup).replace("{{bamf}}", wgse.BAM.disp_FBS)
    result_window_simple(wgse.lang.i18n['MitochondrialDNA'], report, "mHaplo", '14', top=True)

    mainwindow_resume()


def button_yAndMt_BAM():
    """ Processing user button to generate Y and mtDNA only BAM (e.g. for yFull) """
    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN

    # Todo Important to retain same sort order with header for samtools depth processing in the future
    chromos = f'{wgse.BAM.get_chr_name("Y")} {wgse.BAM.get_chr_name("M")}'
    if wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Chr":
        chromos += f' chrY_KI270740v1_random'       # Only Y option; added only to Build 38 based on HG model
    out_qFN = f'"{wgse.outdir.FPB}_chrYM.bam"'

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {out_qFN}\n'

    run_bash_script("ButtonYandMT", commands)

    if not os.path.isfile(nativeOS(unquote(out_qFN))):
        wgse_message("error", 'errYTitle', False, 'errYandMtFile')
    mainwindow_resume()


def button_yOnly_BAM(internal=False):
    """
    Processing user button to generate Y only BAM (always BAM, not CRAM)
    Put in output directory; does not replace currently selected BAM (extract button; not BAM settings button)
    Can be called internally. If so, put in temp area.
    """
    import shutil

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else  ""
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return False     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN
    outybam_qFN = f'"{wgse.outdir.FPB}_chrY.bam"' if not internal else \
                  f'"{wgse.tempf.FP}{wgse.BAM.file_FB}_chrY.bam"'

    outybam_oFN = nativeOS(unquote(outybam_qFN))
    if wgse.BAM.Yonly and not wgse.BAM.Monly:
        # If already a true Yonly BAM, simply copy it.  Internal call requires new name.
        shutil.copy2(wgse.BAM.file_oFN, outybam_oFN)
        return True

    chromos = wgse.BAM.get_chr_name('Y')
    if wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Chr":
        chromos += f' chrY_KI270740v1_random'   # Only Y alt contig in Build 38 model; everyone seems to use it so ...

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {outybam_qFN} \n'

    run_bash_script("ButtonYonly", commands)

    mainwindow_resume() if not internal else ""
    if not os.path.isfile(outybam_oFN):
        wgse_message("error", 'errYTitle', False, 'errYBAMFile')
        return False

    return True


def button_yOnly_VCF():
    """
    Special button to create subset of VCF for Y only; used in Cladefinder.yseq.net.
    We keep the final VCF uncompressed as that is what Cladefinder can deal with.
    """

    bcftools = wgse.bcftoolsx_qFN
    tabix = wgse.tabixx_qFN
    cpus = wgse.os_threads

    # Get correct Reference VCF file based on build, chromosome names and type (All, Yonly)
    refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, type="Yonly")
    if "error" in refVCFtab_qFN or not os.path.isfile(nativeOS(unquote(refVCFtab_qFN))):
        wgse_message("error", 'errYTitle', False, 'errYNoVCFRef')
        mainwindow_resume()
        return

    refgenome_qFN = wgse.BAM.Refgenome_qFN
    if wgse.reflib.missing_refgenome(refgenome_qFN):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist

    # Subset to Y-only BAM if not already available; big time saver when creating called variances
    bamfile_FN = f'{wgse.outdir.FP}{wgse.BAM.file_FB}_chrY.bam'
    bamfile_oFN = nativeOS(bamfile_FN)
    if wgse.BAM.Yonly:
        bamfile_FN = wgse.BAM.file_FN   # Already a subsetted BAM; just use directly
    elif not (os.path.exists(bamfile_oFN) and os.path.getmtime(bamfile_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        # Have to create a sub-setted y-only BAM; use Temp since buried in this routine
        bamfile_FN = f'{wgse.tempf.FP}{wgse.BAM.file_FB}_chrY.bam'   # As created by button_yOnly() in internal mode
        bamfile_oFN = nativeOS(bamfile_FN)
        if not (os.path.exists(bamfile_oFN) and os.path.getmtime(bamfile_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
            # Just check to see if already in temp directory first; otherwise create
            button_yOnly_BAM(internal=True)       # New BAM is in temp due to internal mode call

    # So use bamfile name chosen (without extension) to create various intermediate temp and final file names
    bamfile_FB = os.path.basename(os.path.splitext(bamfile_FN)[0])
    called_bcf = f'"{wgse.tempf.FP}{bamfile_FB}.called.bcf.gz"'             # Note: in temp area and not compressed
    annotated_vcf = f'"{wgse.outdir.FP}{bamfile_FB}.annotated.vcf.gz"'
    bamfile_qFN = f'"{bamfile_FN}"'

    ploidy = "1"    # Y only; not enabled for female so ...
    anncol = "+ID,+INFO/HG,+INFO/ISOGG"

    """  Trying to include and left-justify (normalize) InDels but did not get it working quite right
    temp_norm_bcf = f'"{wgse.tempf.FP}{bamfile_FB}.norm.bcf.gz"'
    f'{bcftools} norm -f {refgenome_qFN} {temp_called_bcf} -Ob -o {temp_norm_bcf}\n'
    f'{bcftools} index {temp_norm_bcf}\n'
    """
    # Note: -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
    # Note: assumes ran "bcftools norm -m +any refVCFTab_qFN" so multiple names merged into single record.
    #        yBrowse has as many as 4 different entries for the same position. Also that CHROM names match.
    # Note: Cladefinder cannot handle BCF so output as compressed VCF
    commands = (
        f'{bcftools} mpileup -B -I -C 50 -T {refVCFtab_qFN} -f {refgenome_qFN} -Ou {bamfile_qFN} | '
        f'  {bcftools} call --ploidy {ploidy} -V indels -m -P 0 --threads {cpus} -Ob -o {called_bcf}\n'
        f'{bcftools} index {called_bcf}\n'
        f'{bcftools} annotate -a {refVCFtab_qFN} -c {anncol} {called_bcf} --threads {cpus} -Oz -o {annotated_vcf}\n'
        f'{tabix} -p vcf {annotated_vcf}\n'
    )
    run_bash_script("AnnotatedVCF-yOnly", commands)

    if not os.path.exists(nativeOS(unquote(annotated_vcf))):
        wgse_message("error", 'errYTitle', False, 'errYVCFFile')
    mainwindow_resume()


def button_annotate_BigYVCF():
    """
    Special for handling FTDNA BigY VCF which comes uncompressed and with incorrect format.
    Must filter first to correct format before using bcftools annotation.
    """
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "Annotate BigY VCF"))
    '''
    # IN DEVELOPMENT; Pseudo-code so far
    if zip compressed:      # FTDNA delivers a .zip file with the source VCF, a readme.txt and regions.bed file
        f'{unzip} {source_VCF} variants.vcf'
        real_VCF = "variants.vcf"
    elif bgzip'ed:
        f'{bgzip} -d {source_VCF}'
        real_VCF = source_VCF.split('.gz')[0]       # Removed trailing .gz on source_VCF name
    else:
        real_VCF = source_VCF
    f'{python} fixFTDNAvcf.py < {source_VCF} > {fixed_VCF}'
    f'{bgzip} {fixed_VCF} ; {bcftools} index {fixed_VCF.gz}'
    f'{bcftools} annotate -a {refVCFtab_qFN} -c +ID,+ID/HG,+ID/ISOGG {fixed_VCF.gz} --threads {cpus} -Ov -o {annot_VCF}'
    f'{bcftools} index {annot_VCF}'
    if compressed:
        f'{mv} {annot_VCF} variants.vcf'
        f'{zip} -f {source_VCF} variants.vcf'
        f'{rm} -f variants.vcf'
    elif bgzip'ed:
        f'{wgsebgzip} {source_VCF}'
    f'{rm} -f {fixed_VCF}'
    '''


def button_ydna_haplogroup():
    """ Processing user button to run yleaf to generate Y Haplogroup call. """

    outd_FP = f'{wgse.tempf.FP}tempYleaf/'
    outf = f'"{outd_FP}haplogroups.txt"'        # note: _qFN
    outd_oFP = nativeOS(outd_FP)
    if os.path.isdir(outd_oFP):             # REH 10Mar2020 If previous run crashed ...
        # wgse.tempf.list.append(outd_oFP)    # Calling clean assures deletion even if DEBUG_MODE is True
        wgse.tempf.clean_item(outd_oFP, ignore_debug_mode=True)

    position_FBS = "WGS_hg38.txt" if "38" in wgse.BAM.Refgenome else "WGS_hg19.txt"
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return  # Check routine reports error if reference does not exist

    filespec = f'-f {wgse.BAM.Refgenome_qFN} -cram {wgse.BAM.file_qFN}' if wgse.BAM.file_type == "CRAM" else \
               f'-bam {wgse.BAM.file_qFN}'

    # Each file path needs to be surrounded by quotes       # REH 10Mar2020
    ylef = f'"{wgse.yleaf_FP}yleaf.py"'
    posf = f'"{wgse.yleaf_FP}Position_files/{position_FBS}"'
    pyth = wgse.python3x_qFN
    samt = wgse.samtoolsx_qFN
    # Todo REH 10Mar2020 bug from wgsev2?  -r specified twice: -r 3 and -r 1?
    # Modified yleaf to accept python and samtools exec path;

    commands = f'{pyth} {ylef} {filespec} -pos {posf} -out {outf} -r 1 -q 20 -b 90 -py {pyth} -samt {samt}\n'

    run_bash_script('ButtonYHaplo', commands)

    result_oFN = nativeOS(unquote(outf))
    if not os.path.exists(result_oFN):
        print('--- FAILURE in yLeaf call (no result file)!')    # Todo convert to wgse_message
        mainwindow_resume()
        return

    with open(result_oFN) as fp:    # MB Jan2020
        fp.readline()               # REH 13Mar2020 This would work if file named correctly; skip header
        fpdata = fp.readline()      # REH 13Mar2020 Read only first result; only supplying one file
        columns_yleaf_result = fpdata.split("\t")

    if len(columns_yleaf_result) < 2 or columns_yleaf_result[1] == 'NA':
        print('--- FAILURE in yLeaf call (NA result)!')    # Todo convert to wgse_message
        mainwindow_resume()
        return

    yhg = columns_yleaf_result[1]   # Second column is Y Haplogroup (in YCC Long form)
    yhgs = columns_yleaf_result[2]  # Third column is Y Haplogroup (in YCC Short form; could be a paragroup)

    out_FB  = wgse.BAM.file_FB
    grep    = wgse.grepx_qFN
    awk     = wgse.awkx_qFN
    awkargs = "-F \"\\\"*\\t\\\"*\" '{print $3}'"   # Todo simplify with raw string r""?
    outsuf  = f'{outd_oFP}{out_FB}{wgse.os_slash}{out_FB}.out'
    termsnp = f'{outd_oFP}snps_terminal.txt'
    # Grep through output listing of 64,000+ SNPs and find all with same haplogroup as terminal
    commands = f'{grep} "{yhg}" "{outsuf}" | {awk} {awkargs} > "{termsnp}"\n'

    run_bash_script('ButtonYHaplo2', commands)

    if not os.path.exists(termsnp):
        print('--- Failure in yLeaf post-processing (no terminal SNPs file)!')    # Todo convert to wgse_message
        mainwindow_resume()
        return

    with open(termsnp) as f:
        terminal_snps = [line.rstrip('\n') for line in f]

    result_window_yHaplo(yhg, yhgs, terminal_snps)


def result_window_yHaplo(yhg, yhgs, terminal_snps):
    """ Y Chromosome Haplogroup call results window pop-up. """
    global yHaploResWindow, font

    font = wgse.fonts.table

    left = yhgs.find('-')
    right = yhgs.find('*')

    total_snp_count = len(terminal_snps)
    if total_snp_count > 0:
        full_snp_list = ", ".join(terminal_snps)        # Full SNP list as comma separated string
        snps_textline = ", ".join(terminal_snps[0:total_snp_count if total_snp_count < 3 else 3])
        if total_snp_count > 3:
            snps_textline += " ..."
        example_snp = terminal_snps[0]

    # Must be a paragroup with no SNPs (or error in processing?)
    else:
        snps_textline = yhgs      # Must be a paragroup; so report its parent haplogroup as SNP list
        full_snp_list = ""        # The full SNPs list (if greater than 3 only so empty here)
        example_snp = yhgs[left+1:right] if -1 < left < right else yhgs # Try to extract haplogroup name (SNP) of parent

    majhap = yhg[0]
    yhg += "*" if right > -1 else ""    # Add paragroup designation to long form

    snps_label  = wgse.lang.i18n['SNPsForThisHG'].replace('{{yhg}}', yhg).replace('{{snps}}', snps_textline)
    yFull_label = wgse.lang.i18n['findHgYFullDescription'].replace('{{BspSNP}}', example_snp)
    ftdna_label = wgse.lang.i18n['findOnFTDNADescLabel'].replace('{{BspSNP}}', example_snp).\
        replace('{{url}}', f'{wgse.ftdna_pubytree_url}{majhap}')

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    yHaploResWindow = Toplevel(wgse.window)
    yHaploResWindow.transient()
    yHaploResWindow.title(wgse.lang.i18n['YChromHaplogroup'])
    yHaploResWindow.geometry(wgse.yHgResult_winsize + offset(wgse.window.geometry(), 50))
    yHaploResWindow.maxsize(width=wgse.yHgResult_maxw, height=wgse.yHgResult_maxh)
    yHaploResWindow.protocol("WM_DELETE_WINDOW", lambda win=yHaploResWindow: button_close(win))
    yHaploResWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''

    yHaploResWindow.rowconfigure(1, weight=1)
    yHaploResWindow.columnconfigure(0, weight=1)
    # On Mac and Linux, window is not sized yet
    wraplength1 = (wgse.yHgResult_maxw - 20) if yHaploResWindow.winfo_width() < 300 else \
        (yHaploResWindow.winfo_width() - 20)
    wraplength2 = int(wraplength1 * 2 / 3)

    # Upper left Haplogroup box
    yhgresFrame = LabelFrame(yHaploResWindow, text=wgse.lang.i18n['FrameHgISOGG'], font=font['16b'])
    yhgresFrame.grid(row=0, column=0, padx=10, pady=5, sticky=NSEW)

    Label(yhgresFrame, text=f"File: {wgse.BAM.disp_FBS}", justify=LEFT,
          font=font['14']).grid(row=0, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=wgse.lang.i18n['YDNAHaploReport'], justify=CENTER, wraplength=wraplength2,
          font=font['14']).grid(row=1, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=f"{yhg}, {yhgs}", justify=CENTER,
          font=font['14b']).grid(row=2, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=font['12']).grid(row=3, column=0, columnspan=4, padx=0, pady=0)

    # Upper right SNPs box
    snpsForHgFrame = LabelFrame(yHaploResWindow, text=wgse.lang.i18n['FrameHgSNPs'], font=font['16b'])
    snpsForHgFrame.grid(row=0, column=1, padx=10, pady=5, sticky=NSEW)

    Label(snpsForHgFrame, text=snps_label, justify=LEFT, wraplength=snpsForHgFrame.winfo_width()-10,
          font=font['14']).grid(row=0, column=0, padx=5, pady=2)
    if total_snp_count > 3:
        Button(snpsForHgFrame, text=wgse.lang.i18n['showAllSnps'], justify=CENTER,
               command=lambda ys=yhg + " SNPs", fs=full_snp_list: result_window_simple(ys, fs, "yHaplo_list", '14'),
               font=font['14']).grid(row=1, column=0, padx=5, pady=2)
    Label(snpsForHgFrame, text=wgse.lang.i18n['SNPsForThisHG2'], justify=CENTER,
          font=font['12']).grid(row=2, column=0, padx=5, pady=2)

    # Lower, spanning-columns Trees box
    findHgOtherTreesFrame = LabelFrame(yHaploResWindow, text=wgse.lang.i18n['FrameHgOtherTrees'], font=font['16b'])
    findHgOtherTreesFrame.grid(row=1, column=0, columnspan=2, padx=10, pady=5, sticky=NSEW)

    DEBUG(f"Wraplengths: window:{wraplength1}, 2/3 frame:{wraplength2}")
    crow = 0
    Label(findHgOtherTreesFrame, text=wgse.lang.i18n['findHgISOGGDescription'], justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnISOGGButtonText'], justify=CENTER,
           command=lambda: webbrowser.open_new(wgse.isogg_tree_url),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Label(findHgOtherTreesFrame, text=yFull_label, justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnYFullButtonText'], justify=CENTER,
           command=lambda: webbrowser.open_new(wgse.yfull_searchsnp_url),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow +=1
    Label(findHgOtherTreesFrame, text=ftdna_label, justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnFTDNAButtonText'], justify=CENTER,
           command=lambda ur=f'{wgse.ftdna_pubytree_url}{yhg[0]}': webbrowser.open_new(ur),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Label(findHgOtherTreesFrame, text=wgse.lang.i18n['warningOtherTreesAreDeeper'], justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1

    Button(yHaploResWindow, text=wgse.lang.i18n['SaveWindow'],
           command=lambda win=yHaploResWindow: button_save(win, "yHaplo"), justify=CENTER,
           font=font['14']).grid(row=2, column=0, padx=5, pady=2)
    Button(yHaploResWindow, text=wgse.lang.i18n['CloseWindow'],
           command=lambda win=yHaploResWindow: button_close(win), justify=CENTER,
           font=font['14']).grid(row=2, column=1, padx=5, pady=2)

    yHaploResWindow.update()
    yHaploResWindow.grab_set()      # Toplevel equivalent of mainloop
    yHaploResWindow.wait_window()


def result_window_simple(title, simple_text, file_ext, pts, top=False, wwrap="word"):
    """
    Generic simple result window handler.  Provides results given in simple_text and a save and exit button.
    """
    global font

    font = wgse.fonts.table

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    simResWindow = Toplevel(wgse.window)
    simResWindow.transient()
    simResWindow.title(title)
    simResWindow.geometry(wgse.simResult_winsize + offset(wgse.window.geometry(), 50))
    simResWindow.protocol("WM_DELETE_WINDOW", lambda win=simResWindow: button_close(win))
    simResWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    simResWindow.maxsize(width=wgse.simResult_maxw, height=wgse.simResult_maxh)

    simResWindow.rowconfigure(1, weight=1)
    simResWindow.columnconfigure(0, weight=1)
    # On Mac and Linux, window is not sized yet
    # wraplength1 = (wgse.simResult_maxw - 20) if simResWindow.winfo_width() < 300 else \
    #     (simResWindow.winfo_width() - 20)

    Label(simResWindow, text=f'{title}:', justify=LEFT,
          font=font['14b']).grid(row=0, column=0, columnspan=3, padx=5, pady=2)

    # Scrollable text body with main content
    bodyFrame = Frame(simResWindow)
    bodyFrame.grid(row=1, column=0, columnspan=3, padx=1)
    bodyFrame.rowconfigure(0, weight=1)
    bodyFrame.columnconfigure(0, weight=1)

    bodyText = Text(bodyFrame, width=wgse.simResult_maxw, height=wgse.simResult_maxh, wrap=wwrap, font=font[pts])
    bodyText.insert('end', simple_text) # justify=LEFT, wraplength=wraplength1,
    bodyText.grid(row=0, column=0, padx=5, pady=2, sticky=NSEW)

    bodyScrolly = Scrollbar(bodyFrame, orient=VERTICAL, command=bodyText.yview)
    bodyScrolly.grid(row=0, column=1, padx=2, pady=2, sticky=NS)
    bodyText.config(yscrollcommand=bodyScrolly.set)

    if wwrap == "none":      # Need X scrollable also
        bodyScrollx = Scrollbar(bodyFrame, orient=HORIZONTAL, command=bodyText.xview)
        bodyScrollx.grid(row=1, column=0, padx=2, pady=2, sticky=EW)
        bodyText.config(xscrollcommand=bodyScrollx.set)

    # Bottom row information and buttons
    Label(simResWindow, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=font['12']).grid(row=2, column=0, columnspan=3, padx=0, pady=0)

    Button(simResWindow, text=wgse.lang.i18n['SaveWindow'],
           command=lambda win=simResWindow, fe=file_ext: button_save(win, fe), justify=CENTER,
           font=font['14']).grid(row=3, column=0, padx=5, pady=2)
    Button(simResWindow, text=wgse.lang.i18n['CloseWindow'],
           command=lambda win=simResWindow, t=top: button_close(win, t), justify=CENTER,
           font=font['14']).grid(row=3, column=1, padx=5, pady=2)
    Label(simResWindow, text="   ", font=font['12']).grid(row=3, column=2, padx=5, pady=2)

    simResWindow.update()
    simResWindow.grab_set()     # Toplevel equivalent of mainloop
    simResWindow.wait_window()


def button_export_unmapped_reads(internal=False):
    """ Process user button to export unmapped reads from BAM. """
    # Todo include key non-human decoys that should be considered unmapped (EBV, hs37d5?, etc)

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist

    paired = wgse.BAM.ReadType == "Paired"

    samtools = wgse.samtoolsx_qFN
    bamfile = wgse.BAM.file_qFN
    tempdir = f'"{wgse.tempf.FP}"'
    commands = ""

    # Subset BAM to just unmapped; if not already existing
    unmap_qFN = f'"{wgse.outdir.FPB}_unmap.bam"'
    unmap_oFN = f'{wgse.outdir.oFPB}_unmap.bam'
    unmap_good = os.path.exists(unmap_oFN) and os.path.getsize(unmap_oFN) > 10000 and \
                 os.path.getmtime(unmap_oFN) > os.path.getmtime(wgse.BAM.file_oFN)

    if not unmap_good:
        commands += (
            f"{samtools} view -bh {cram_opt} -@ {wgse.os_threads} -o {unmap_qFN} {bamfile} '*' \n"
            f"{samtools} index {bamfile} \n"
        )

    # Sort unmapped BAM in name order, then call fastq command to (split and) write FastQ's
    r1_qFN = r2_qFN = r_qFN = None
    files_exist = False
    if paired:
        r1_qFN = f'"{wgse.outdir.FPB}_unmapR1.fastq.gz"'
        r2_qFN = f'"{wgse.outdir.FPB}_unmapR2.fastq.gz"'
        r1_oFN = f'{wgse.outdir.oFPB}_unmapR1.fastq.gz'
        r2_oFN = f'{wgse.outdir.oFPB}_unmapR2.fastq.gz'
        if unmap_good:
            r1_good = os.path.exists(r1_oFN) and os.path.getmtime(r1_oFN) > os.path.getmtime(unmap_oFN)
            r2_good = os.path.exists(r2_oFN) and os.path.getmtime(r2_oFN) > os.path.getmtime(unmap_oFN)
            files_exist = unmap_good and r1_good and r2_good
    else:
        r_qFN = f'"{wgse.outdir.FPB}_unmap.fastq.gz"'
        r_oFN = f'{wgse.outdir.oFPB}_unmap.fastq.gz'
        if unmap_good:
            r_good = os.path.exists(r_oFN) and os.path.getmtime(r_oFN) > os.path.getmtime(unmap_oFN)
            files_exist = unmap_good and r_good

    if not files_exist:
        # Unmapped should normally be very small.  But just in case it is larger, pretend the unmapped file is 33%
        #  the size of the BAM. Need to make sure enough CPUs / memory to do the name sort on MacOS
        sort_mem, sort_cpus = _adjust_mem_threads(wgse.BAM.file_stats.st_size // 3, wgse.BAM.file_type, "Name")
        if sort_cpus == 0:  # Not enough memory to run samtools sort on MacOS; already reported in _adjust call
            mainwindow_resume()
            return False

        # bam2fq deprecated in 1.9; now use fastq
        files = f"-1 {r1_qFN} -2 {r2_qFN} -0 /dev/null -s /dev/null" if paired else \
                f"-1 /dev/null -2 /dev/null -0 {r_qFN} -s {r_qFN}"
        commands += (
            f"{samtools} sort -n -T {tempdir} -m {sort_mem} -@ {sort_cpus} {unmap_qFN} | "
            f"  {samtools} fastq {files} -n -@ {wgse.os_threads}\n"
        )

    if commands:
        run_bash_script('ButtonUnmappedReads', commands)

    if internal:        # internal call; not from a user button click. So delete the intermediary unmap BAM file
        wgse.tempf.list.append(unmap_oFN)
    else:               # User button click, so give a result window pop-up
        # Want just file names; not full path
        unmapbam = f'{wgse.BAM.file_FB}_unmap.bam'
        script = wgse.lang.i18n['DescriptionUploadCosmosId'].replace("{{obam}}", unmapbam)

        if paired:
            unmapR1fastq = f'{wgse.BAM.file_FB}_unmapR1.fastq.gz'
            unmapR2fastq = f'{wgse.BAM.file_FB}_unmapR2.fastq.gz'
            tscript = script.replace("{{fq1}}", unmapR1fastq).replace("{{fq2}}", unmapR2fastq)
        else:
            unmapfastq = f'{wgse.BAM.file_FB}_unmap.fastq.gz"'
            tscript = script.replace("{{fq1}}", unmapfastq).replace("{{fq2}}", "(single-end)")

        result_window_simple(wgse.lang.i18n['HowToContinue'], tscript, "unmap", '14', top=True)

    mainwindow_resume()


def button_extract_WES():
    """
        Really only needed for those who did a WGZ with Dante (30x WGS and 130x WES). But helpful to get a small,
        more workable file of Exome-only for processing to VCF and studying known SNPs for health.
        If a Y or Y and MT only; use special BED file not based on Exomes but Y CombBED / McDonald / Poznik region
    """

    outfile_qFN = f'"{wgse.outdir.FPB}_wes{wgse.BAM.file_FS}"'
    outfile_oFN = nativeOS(unquote(outfile_qFN))

    samtools = wgse.samtoolsx_qFN
    bamfile = wgse.BAM.file_qFN
    bedfile = wgse.reflib.get_wes_bed_file_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, wgse.BAM.Yonly)
    cramopts = f'-C --reference {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else "-b"

    if (bedfile and not check_exists(nativeOS(unquote(bedfile)), 'errNoBEDFile')) or \
       wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return

    commands = (
        f'{samtools} view -h -L {bedfile} {cramopts} -@ {wgse.os_threads} -o {outfile_qFN} {bamfile} \n'
        f'{samtools} index {outfile_qFN} \n'
    )

    run_bash_script("ExtractWES", commands, parent=wgse.window)

    # Check if WES BAM file generated
    if not os.path.isfile(outfile_oFN):
        # Todo report error generating via wgse_message before returning
        mainwindow_resume()
        return

    set_BAM_file(unquote(outfile_qFN))    # Replace selected BAM and garbage collect old one
    mainwindow_resume()


minVCFSize = 50000


def button_SNP_VCF(gVCF=False):
    _button_VCF(vtype="SNP", gVCF=gVCF)


def button_InDel_VCF(gVCF=False):
    _button_VCF(vtype="InDel", gVCF=gVCF)


def _button_VCF(vtype="SNP", gVCF=False):

    bcftools = wgse.bcftoolsx_qFN
    tabix = wgse.tabixx_qFN
    cpus = wgse.os_threads
    ploidy_qFN = f'"{wgse.reflib.cma_FP}ploidy.txt"'

    gvcf_opt = f"-g 1" if gVCF else "-v"  #  add -M to keep masked N sites? add -i N to put reference value if pileup less than N?
    vtype = vtype.lower()
    notv = "-V indels" if vtype == "snp" else \
           "-V snps" if vtype == "indel" else \
           "" if vtype == "both" else ""

    final_vcf_FN  = f'{wgse.outdir.FPB}.{vtype}.vcf.gz'
    final_vcf_oFN = nativeOS(final_vcf_FN)
    final_vcf_qFN = f'"{final_vcf_FN}"'
    bamfile = f'"{wgse.BAM.file_FN}"'

    refgen_qFN = wgse.BAM.Refgenome_qFN
    if wgse.reflib.missing_refgenome(refgen_qFN):
        return  # check() reports error if reference genome does not exist

    if not (os.path.exists(final_vcf_oFN) and
            os.path.getmtime(final_vcf_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and
            os.path.getsize(final_vcf_oFN) > minVCFSize):

        commands = (
            f'{bcftools} mpileup {gvcf_opt} -B -I -C 50 -f {refgen_qFN} -Ou {bamfile} | '
            f'  {bcftools} call --ploidy-file {ploidy_qFN} {notv} {gvcf_opt} -m -P 0 --threads {cpus} -Oz -o {final_vcf_qFN}\n'
            f'{tabix} -p vcf {final_vcf_qFN}\n'
        )

        run_bash_script("ButtonVCF", commands, parent=wgse.window)


def button_CNV_VCF():
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "CNV VCF"))


def button_SV_VCF():
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "SV VCF"))


def button_annotate_VCF():
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "Annotated VCF"))


def button_filter_VCF():
    wgse_message("info", 'ComingSoonTitle', True, wgse.lang.i18n['ComingSoonBody'].replace("{{FEATURE}}", "Filtered VCF"))


def _button_unindex_BAM():
    if wgse.BAM and wgse.BAM.Indexed and wgse.BAM.file_type in ["BAM", "CRAM"]:
        bam = wgse.BAM.file_oFN + (".bai" if wgse.BAM.file_type == "BAM" else ".crai")
        wgse.tempf.list.append(bam)
        wgse.BAM.Indexed = False
        set_BAM_window_settings()   # Do not need full new BAM file setup / read-in; just process for Index file gone

    mainwindow_resume()


def _button_unsort_BAM():
    """ Special internal call enabled only in DEBUG_MODE; simply updates the existing BAM Header to unsorted. """
    if wgse.BAM and wgse.BAM.Sorted:
        samtools = wgse.samtoolsx_qFN
        sed = wgse.sedx_qFN
        newhead = f'"{wgse.tempf.FP}newhead.sam"'
        bam = wgse.BAM.file_qFN
        if "_sorted" in bam:    # Try to avoid adding conflicting names
            unsortbam = bam.replace("_sorted", "_unsorted")
        else:                   # Just append unsorted status
            unsortbam = f'"{wgse.BAM.file_FPB}_unsorted.bam"'

        # CRAM decode requires reference genome be specified with VIEW command
        cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
        if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
            mainwindow_resume()
            return  # Check routine reports error if reference does not exist

        # Will only work where sed is in the PATH; cannot use quoted, path-based sed call in Samtools it seems (Win10)
        # cramopts = "-i" if wgse.BAM.file_type == "CRAM" else ""    # will do in place but we want a copy
        # commands = f"""{samtools} reheader {cramopts} -c 'sed "s/coordinate/unsorted/"' {bam} > {unsortbam}\n"""
        commands = (
            f'{samtools} view -H {cram_opt} {bam} | {sed} "s/SO:coordinate/SO:unsorted/" > {newhead} \n'
            f'{samtools} reheader {newhead} {bam} > {unsortbam}\n'
        )

        run_bash_script("UnsortBAM", commands)

        set_BAM_file(unquote(unsortbam))     # Replace current BAM and garbage collect old
    else:
        DEBUG("***Internal ERROR: called unsort_BAM with unsort'ed BAM. How was button available?")

    mainwindow_resume()


def _button_fontsize():
    global fontsizeButton, fontfaceButton, fontsetLabel, wresetButton

    fontsize = simpledialog.askinteger(
        wgse.lang.i18n['FontSizeRequest'], wgse.lang.i18n['FontSizeContent'], parent=wgse.window,
        initialvalue=wgse.fonts.basept, minvalue=0, maxvalue=20)

    if not fontsize or fontsize == 0:
        fontsize = 0
    elif not (6 <= fontsize <= 20):
        wgse_message("error", 'ValueOutOfRangeTitle', True,
                     wgse.lang.i18n['FontSizeOutOfRange'].replace('{{FONTSIZE}}', fontsize))
        return
    elif fontsize == wgse.fonts.basept:
        return          # No change; nothing to do

    DEBUG(f"Base Font Point Size requested: {fontsize}, current: {wgse.fonts.basept}, "
          f"platform default: {wgse.fonts.default_basept} pts")
    wgse.fonts.change(newbasept=fontsize)                   # New user setting; 0 to clear. Changes wgse.fonts.basept

    wgse.save_settings()                                    # Save new setting (clear saved setting if default)

    override = wgse.fonts.default_basept != wgse.fonts.basept
    curfont = (wgse.fonts.face, wgse.fonts.basept+2)        # Nominally use 14 (base 12 + 2)
    fontsizeButton.configure(                               # Update button value and color as necessary
        text=str(wgse.fonts.basept)+' pt', bg=cyan if override else wgse.defbut_bg, font=curfont)

    fontsetLabel.configure(font=curfont)
    fontfaceButton.configure(font=curfont)
    wresetButton.configure(font=curfont)


global _fontfaceWindow, _fontface


def _button_fontface():
    from tkinter.font import families

    global _fontfaceWindow, _fontface

    _fontfaceWindow = Toplevel(wgse.window)
    _fontfaceWindow.transient()
    _fontfaceWindow.title(wgse.lang.i18n["FontFaceRequestTitle"])
    _fontfaceWindow.geometry("450x200" + offset(wgse.window.geometry(), 50))
    _fontfaceWindow.protocol("WM_DELETE_WINDOW", lambda win=_fontfaceWindow: button_close(win, False))
    _fontfaceWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    _fontfaceWindow.columnconfigure(0, weight=1)
    _fontfaceWindow.rowconfigure(1, weight=1)

    Label(_fontfaceWindow, text=wgse.lang.i18n["FontFaceRequest"], font=font['14'], justify="left"
          ).grid(row=0, padx=5, pady=5, sticky=NW)

    _fontface = StringVar(_fontfaceWindow)

    _fontfaceCB = Combobox(_fontfaceWindow, textvariable=_fontface, values=families(), state='readonly',
                           height=10, width=60, justify="left", )
    _fontfaceCB.set(wgse.fonts.face)
    _fontfaceCB.bind('<<ComboboxSelected>>', _button_fontface_handler)
    _fontfaceCB.grid(row=1, sticky=NW)

    _fontfaceWindow.update()
    _fontfaceWindow.grab_set()     # Toplevel equivalent of mainloop
    _fontfaceWindow.wait_window(_fontfaceWindow)


def _button_fontface_handler(event):
    global fontfaceButton, fontsetLabel, fontsizeButton, wresetButton
    global _fontfaceWindow, _fontface

    _fontfaceWindow.destroy()

    newfontface = _fontface.get()

    if newfontface != wgse.fonts.face:
        wgse.fonts.change(newface=newfontface)
        wgse.save_settings()                                    # Save new setting (clear saved setting if default)

        override = wgse.fonts.default_face != wgse.fonts.face
        curfont = (wgse.fonts.face, wgse.fonts.basept+2)        # Nominally use 14 (base 12 + 2)
        flabel = wgse.fonts.face[:22]+'...' if len(wgse.fonts.face) > 25 else wgse.fonts.face
        fontfaceButton.configure(                               # Update button value and color as necessary
            text=flabel, bg=cyan if override else wgse.defbut_bg, font=curfont)

        fontsetLabel.configure(font=curfont)
        fontsizeButton.configure(font=curfont)
        wresetButton.configure(font=curfont)

    mainwindow_resume()


def _button_maxmem_setting():
    global maxmemButton

    maxmem = simpledialog.askinteger(
        wgse.lang.i18n['MaxMemTitle'], wgse.lang.i18n['MaxMemContent'], parent=wgse.window,
        initialvalue=wgse.os_totmem_saved, minvalue=0, maxvalue=100)

    if not maxmem or maxmem == 0:
        maxmem = 0
    elif not (2 <= maxmem < 100):      # Allow any maximum up to 99; will not exceed os_totmem anyway if active
        wgse_message("error", 'ValueOutOfRangeTitle', True,
                     wgse.lang.i18n['MaxMemOutOfRange'].replace('{{maxmem}}', wgse.os_totmem))
        return

    DEBUG(f"Max Tot Mem requested: {maxmem}, saved: {wgse.os_totmem_saved}, "
          f"platform avail: {wgse.os_totmem_proc//10**9}, previously set: {wgse.os_totmem//10**9} GB")

    # If changed user setting, save new value and configure button value
    if maxmem != wgse.os_totmem_saved:
        wgse.os_totmem_saved = maxmem               # New user setting saved (0 means not set)
        wgse.save_settings()                        # Save new settings file (0 clears the setting)
        maxmemButton.configure(text=str(maxmem))    # Change button to new value

    # Depending on whether user override is active or not, calculate appropriate (new) values ...
    override = 2 <= maxmem < wgse.os_totmem_proc // 10**9                   # Is user override active or not
    maxmemButton.configure(bg=cyan if override else wgse.defbut_bg)         # Color button when override
    wgse.os_totmem = maxmem * 10**9 if override else wgse.os_totmem_proc    # Save new value
    wgse.os_mem = wgse.set_mem_per_thread_millions(wgse.os_totmem, wgse.os_threads)     # Note: a string

    DEBUG(f'Updated settings: total Mem: {wgse.os_totmem//10**9} GB, mem per thread: {wgse.os_mem}B')


def _button_maxthread_setting():
    global maxthreadButton

    maxthread = simpledialog.askinteger(
        wgse.lang.i18n['MaxThreadTitle'], wgse.lang.i18n['MaxThreadContent'], parent=wgse.window,
        initialvalue=wgse.os_threads_saved, minvalue=0, maxvalue=100)

    if not maxthread or maxthread == 0:
        maxthread = 0
    elif not (1 <= maxthread < 100):    # Allow any maximum up to 99; will not exceed os_threads anyway if active
        wgse_message("error", 'ValueOutOfRangeTitle', True,
                     wgse.lang.i18n['MaxThreadOutOfRange'].replace('{{maxthread}}', wgse.os_threads))
        return

    DEBUG(f"Max Threads requested: {maxthread}, saved: {wgse.os_threads_saved}, "
          f"platform avail: {wgse.os_threads_proc}, previously set: {wgse.os_threads}")

    # If changed user setting, save new value and configure button value
    if maxthread != wgse.os_threads_saved:
        wgse.os_threads_saved = maxthread                   # New user setting saved (0 means not set)
        wgse.save_settings()                                # Save new settings file (0 clears the setting)
        maxthreadButton.configure(text=str(maxthread))      # Change button to new value

    # Depending on whether user override is active or not, calculate appropriate (new) values ...
    override = 1 <= maxthread < wgse.os_threads_proc                            # If user override active or not
    maxthreadButton.configure(bg=cyan if override else wgse.defbut_bg)          # Color button when active
    wgse.os_threads = maxthread if override else wgse.os_threads_proc           # Save new value
    wgse.os_mem = wgse.set_mem_per_thread_millions(wgse.os_totmem, wgse.os_threads)     # Note: a string

    DEBUG(f'Updated settings: Threads: {wgse.os_threads}, mem per thread: {wgse.os_mem}B')


def _button_SubsetBAM():
    """
        Create BAM subset. Based on percentage parameter, subset a BAM file to make it smaller. Use random subset
        capability with samtools view -s . Auto handled CRAM to CRAM or BAM to BAM.
    """
    subset = simpledialog.askfloat(wgse.lang.i18n['SubsetBAMAskRangeTitle'], wgse.lang.i18n['SubsetBAMAskRangeContent'],
                                   parent=wgse.window)
    DEBUG(f"Subset: {subset}")
    if not (subset and 0.0 <= subset <= 100.0):
        wgse_message("error", 'ValueOutOfRangeTitle', False, 'SubsetOutOfRange')
        return

    samtools = wgse.samtoolsx_qFN

    # CRAM decode requires reference genome be specified with VIEW command
    if wgse.reflib.missing_refgenome(wgse.BAM.Refgenome_qFN, required=wgse.BAM.file_type == "CRAM"):
        mainwindow_resume()
        return      # Check routine reports error if reference does not exist
    ext, cram_opt = ('.cram', f'-C -T {wgse.BAM.Refgenome_qFN}') if wgse.BAM.file_type == "CRAM" else \
                    ('.bam', "-b")


    BAM_qFN = wgse.BAM.file_qFN
    SUB_FN = f'{wgse.outdir.FPB}_{int(subset)}%{ext}'   # Percent in name seems legal everywhere; albeit rare in use
    SUB_qFN = f'"{SUB_FN}"'
    DEBUG(f"Output Subset BAM (final): {SUB_FN}")

    subset /= 100.0
    commands  = (
        f'{samtools} view -bh -s {subset} {cram_opt} -@ {wgse.os_threads} -o {SUB_qFN} {BAM_qFN} \n'
        f'{samtools} index {SUB_qFN} \n'
    )

    run_bash_script("GenBAMSubset", commands)

    set_BAM_file(SUB_FN)    # Replace selected BAM and garbage collect old one

    mainwindow_resume()


global _track, _changeTrackWindow


def _button_change_release_track():
    global _track, _changeTrackWindow

    _changeTrackWindow = Toplevel(wgse.window)
    _changeTrackWindow.transient()
    _changeTrackWindow.title(wgse.lang.i18n["TrackChangeRequestTitle"])
    _changeTrackWindow.geometry("450x200" + offset(wgse.window.geometry(), 50))
    _changeTrackWindow.protocol("WM_DELETE_WINDOW", lambda win=_changeTrackWindow: button_close(win, False))
    _changeTrackWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    _changeTrackWindow.columnconfigure(0, weight=1)
    _changeTrackWindow.rowconfigure(1, weight=1)

    _tracks = [wgse.lang.i18n["Beta"], wgse.lang.i18n["Alpha"], wgse.lang.i18n["Dev"]]

    Label(_changeTrackWindow, text=wgse.lang.i18n["TrackChangeRequest"], font=font['14'], justify="left"
          ).grid(row=0, padx=5, pady=5, sticky=NW)

    _track = StringVar(_changeTrackWindow)

    _changeTrackCB = Combobox(_changeTrackWindow, textvariable=_track, values=_tracks, state='readonly',
                              height=10, width=60, justify="left")
    _changeTrackCB.set(wgse.lang.i18n[wgse.track])
    _changeTrackCB.bind('<<ComboboxSelected>>', _button_changeTrack_handler)
    _changeTrackCB.grid(row=1, sticky=NW)

    _changeTrackWindow.update()
    _changeTrackWindow.grab_set()  # Toplevel equivalent of mainloop
    _changeTrackWindow.wait_window(_changeTrackWindow)


def _button_changeTrack_handler(event):
    global releaseTrackButton, _track, _changeTrackWindow

    _changeTrackWindow.destroy()

    # Rertieve user selection and convert back into internal form (from i18n)
    i18n_track = _track.get()
    newtrack = "Beta"  if i18n_track == wgse.lang.i18n["Beta"] else \
               "Alpha" if i18n_track == wgse.lang.i18n["Alpha"] else \
               "Dev"   if i18n_track == wgse.lang.i18n["Dev"] else None

    if newtrack is not None and newtrack != wgse.track:
        update_release_track(newtrack)
        load_wgse_version()  # Create new program version string; updates wgse.track also
        releaseTrackButton.configure(text=i18n_track)

        DEBUG(f'WGS Extract: {wgse.__version__}')
        version_fg = 'red' if is_wgse_outdated() else wgse.deflab_fg
        versionLabel.configure(text=wgse.__version__, fg=version_fg)

    mainwindow_resume()


def button_debug_mode_toggle():
    """ Toggle debug mode.  Note, as not in the settings file, must add or remove the file for startup mode. """
    from pathlib import Path

    if wgse.DEBUG_MODE:
        DEBUG("***** DEBUG_MODE off *****")
        os.remove(wgse.debugset_oFN)
        wgse.DEBUG_MODE = False
    else:
        wgse.DEBUG_MODE = True
        Path(wgse.debugset_oFN).touch()
        DEBUG("***** DEBUG_MODE on  *****")

    mainwindow_reset()
    mainwindow_resume()


def button_language_reload():
    """
    Setup in DEBUG_MODE section for developers to reload the languages.csv file (presumably after changes)
    and recreate the mainwindow like occurs after a language change. This along with language setting button allows
    for quick loops of language translation development (test by quick change and reload)
    """
    from utilities import LanguageStrings       # Localized import to prevent module import loop

    wgse.save_settings()        # Make sure we have saved the latest settings; like on exit

    # Restart the language subsystem by doing a new class init and assigning to global setting; garbage collect old
    wgse.lang = LanguageStrings(wgse.language_oFN)  # (re)Start language subsystem (utilities.py); reread language.csv

    wgse.load_settings(False)   # re-setup selected language; do not reload BAM (not top level)

    mainwindow_resume()


def button_change_language():
    """
    User button that displays currently set language. But when clicked, requests the new language setting from the user
    It then proceeds to do a switch_language setting which resets the language setting class and causes a new
    mainWindow to be created and drawn for the new language
    """
    # Todo move GUI / get_language from utilities/Language to this module. Then simply call get_language.
    #   Call get_language in this module if settings does not include language at startup like before.
    button_set_language()     # Do everything locally in the language class

    mainwindow_resume()


def button_set_language():
    """
       Language subsystem user pop-up query for language setting to use. Note, if no stored language in settings,
       then this is called before anything else so we have a language for the initial program window. Hence always
       creating our own root window.
    """
    global selectLanguageWindow, font

    # Font subsystem may not be setup yet; so have reasonable default
    font = wgse.fonts.table if wgse and wgse.fonts else {'14': ("Times New Roman", 14)}

    DEBUG("Displaying Language Selection Pop-up")
    if wgse.window:             # Withdraw main window as will destroy and create new one with new language
        wgse.window.withdraw()
        selectLanguageWindow = Toplevel(wgse.window)
        selectLanguageWindow.transient()
        selectLanguageWindow.geometry(offset(wgse.window.geometry(), 50))
    else:
        selectLanguageWindow = Tk()
        selectLanguageWindow.geometry("")

    selectLanguageWindow.title(wgse.lang.i18n["RequestLanguage"] if wgse.lang.language else "Please Select Language")
    selectLanguageWindow.protocol("WM_DELETE_WINDOW", lambda win=selectLanguageWindow: button_close(win, False))
    selectLanguageWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    selectLanguageWindow.columnconfigure(0, weight=1)
    selectLanguageWindow.rowconfigure(0, weight=1)

    crow = 0
    for lang in wgse.lang.avail_langs:       # Need all languages so go back to original class dictionary
        Label(selectLanguageWindow, text=wgse.lang.request[lang],
              font=font['14']).grid(column=0, row=crow, padx=5, pady=2)
        Button(selectLanguageWindow, text=lang, font=font['14'],
               command=lambda l=lang: wgse.lang.switch_language(l)).grid(column=1, row=crow, padx=5, pady=2)
        crow += 1

    wgse.lang.selectLanguageWindow = selectLanguageWindow   # So we can close window in switch_language() after button
    selectLanguageWindow.update()
    selectLanguageWindow.grab_set()     # Toplevel equivalent of mainloop
    if wgse.window:           # If wgse.window set, have run mainwindow_setup
        selectLanguageWindow.wait_window(selectLanguageWindow)
        # mainwindow_resume()       # Handled in switch_language()
    else:
        selectLanguageWindow.mainloop()  # Run because we ask for language before setting up main window and its loop


def button_set_referencelibrary():
    """
    Button to allow override of default Reference_Library location.  Note: user must move content or rerun Library
    command after setting to reinstall the reference library.
    """
    global reflibDirectoryButton

    initialdir = os.path.dirname(wgse.reflib.FP[:-1]) if wgse.reflib and wgse.reflib.FP else \
                 os.path.dirname(wgse.reflib.default_FP[:-1]) if wgse.reflib and wgse.reflib.default_FP else \
                 wgse.install_FP
    reflib_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectReferenceLibrary'],
                                        initialdir=initialdir)  # Returns Unix/Universal, and no trailing slash
    if not reflib_FP:      # If still not set ...
        DEBUG(f"New Reference Library Directory not set (cancelled)")
    else:
        reflib_FP += '/' if reflib_FP[-1] != '/' else ''    # Note: not NativeOS so forward slash always

        wgse.reflib.change(reflib_FP)
        reflibDirectoryButton.configure(text=wgse.reflib.FB if wgse.reflib.set else wgse.lang.i18n['Default'])
        reflibDirectoryButton.configure(font=font['14b'] if wgse.reflib.set else font['14'])
        reflibDirectoryButton.configure(bg=cyan if wgse.reflib.set else wgse.defbut_bg)

    mainwindow_resume()


def button_set_tempdir():
    """
    Button to allow override of default Temporary Files location.  Note: user selects and must already exist.
    No files to check for as should be empty to start.  Will clean out old one before replacing.
    """
    global tempDirectoryButton

    initialdir = os.path.dirname(wgse.tempf.FP[:-1]) if wgse.tempf and wgse.tempf.FP else \
                 os.path.dirname(wgse.tempf.default_FP[:-1]) if wgse.tempf and wgse.tempf.default_FP else \
                 wgse.install_FP
    temp_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectTemporaryFiles'],
                                      initialdir=initialdir)  # Returns Unix/Universal, and no trailing slash
    if not temp_FP:  # If still not set ...
        # User likely hit exit or cancel; simply return (do not clear old value)
        DEBUG(f"New Temporary File Directory not set (cancelled)")
    else:
        # Assure trailing slash
        temp_FP += '/' if temp_FP[-1] != '/' else ''    # Note: universalOS so forward slash always

        wgse.tempf.clean(True)                      # Clean old area before leaving
        wgse.tempf.change(temp_FP)                  # Setup new Temporary Files area (checking validity first)
        tempDirectoryButton.configure(text=wgse.tempf.FB if wgse.tempf and wgse.tempf.set else
                                      wgse.lang.i18n['Default'])

    mainwindow_resume()


def button_toggle_preferred_server():
    """ Toggle current value of the Preferred server for automatic reference genome downloads."""
    global prefServerButton

    wgse.prefserver = "EBI" if wgse.prefserver == "NIH" else "NIH"
    prefServerButton.config(text=wgse.prefserver)
    mainwindow_resume()


def button_wsl_bwa_patch():
    """ Toggle wsl_bwa_patch setting on or off.  True is On / Active / use WSL BWA."""
    global wslbwaButton

    wgse.wsl_bwa_patch = not wgse.wsl_bwa_patch
    wslbwaButton.configure(text=wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else wgse.lang.i18n['Inactive'])
    mainwindow_resume()


def mainwindow_init():
    """
    mainWindow (all GUI processing) subsystem.  Not yet a class but preparing for that.  So this is the early __init__.
    withdraw until cal mainwindow_setup to fill it in later
    """
    global font

    # Create root / main window so toplevel dialogs can be made
    root = Tk()         # Like for mainloop, should only be one in whole application. This is it.
    root.title("WGS Extract")
    root.geometry(wgse.mainwindow_size)
    root.protocol("WM_DELETE_WINDOW", button_exit)
    root.columnconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)

    # Let's hide window until mainwindow_setup fills it in. Created now so can have dialogs like language selection.
    root.withdraw()

    # iconbitmap call causes window to appear. If before withdraw, causes window to flash on
    root.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''

    return root


def mainwindow_reset():
    if wgse and wgse.window and wgse.dnaImage:  # If main window setup, then recreate with new language and DEBUG_MODE
        DEBUG('Resetting / redrawing main window')
        wgse.dnaImage = None  # Only way we know mainwindow_setup was run
        wgse.window.destroy()  # Kill current mainwindow
        top_level = False
        wgse.window = mainwindow_init()         # Create new mainWindow in new language or fonts settings
        wgse.fonts = FontTypes(top_level)       # Create new Fonts subsystem (requires window setup first)
        wgse.load_settings(top_level)           # Reload settings as may override default fonts, etc.
        mainwindow_setup()                      # Fill in new mainWindow (sets wgse.dnaImage)
        # Need to call mainloop again after window destroy?


def mainwindow_setup():
    """
    Setup of TK main window with all the pane's, buttons and labels. Will be looped over for rest of program life.
    Could be renamed to change_mainWindow() to be similar to other classes.  Only called again when recreating
    after a language change.

    wgse.window not null says mainwindow_init() call made.  wgse.dnaImage not null says this mainwindow_setup called.
    """
    global all_action_buttons, font

    # Need to globally reset labels and buttons for Settings (Output Dir, BAM File)
    global versionLabel, documentationButton, exitButton, titleFileLabel, titleCommandLabel
    global outdirLabel, outdirButton, outdirUnselectButton, langSelectLabel, langSelectButton, lreloadButton
    global reflibDirectoryLabel, reflibDirectoryButton, tempDirectoryLabel, tempDirectoryButton
    global prefServerButton, bamSelectedLabel, bamSelectButton, bamReferenceGenomeLabel
    # Need to globally adjust summary stats of BAM (and enable/disable its detailed stats button)
    global infoLabelMapAvgReadDepth, infoLabelGender, infoLabelChroms
    global bamMapAvgReadDepthLabel, bamGenderLabel, bamChromsLabel, bamFsizeLabel      # bamAverageReadLengthLabel,
    global bamIdxstatsButton, bamHeaderButton, bamIndexButton, bamSortButton, bamConvertButton, bamWESButton
    global bamRealignButton, bamUnselectButton
    # Need to globally adjust the action buttons on Extract and Analysis tabs
    global autosomalFormatsButton
    global mitoFASTAButton, mitoBAMButton, mitoVCFButton
    global yANDmtButton, yOnlyButton, yVCFButton
    global haplogroupYButton, haplogroupMtButton, exportUnmappedReadsButton
    global bamAlignButton, bamUnalignButton, fastqFastpButton, fastqFastqcButton
    global SNPVCFButton, InDelVCFButton, CNVVCFButton, SVVCFButton, AnnotateVCFButton, FilterVCFButton, VarQCButton
    # Debug_MODE only buttons
    global wslbwaButton, runMicroParallelButton, subsetBAMButton
    global maxsettingsLabel, maxmemButton, maxmemLabel, maxthreadButton, maxthreadLabel
    global fontsetLabel, fontsizeButton, fontfaceButton, wresetButton, releaseTrackButton

    def_bg = wgse.window.cget("background")

    # Only for MacOS Python Tkinter; they alternate default background colors to nake all mested widgets visible
    wgse.style = Style()
    wgse.style.configure('TFrame', background=def_bg)    # '#F0F0F0' 'SystemButtonFace'  '#d9d9d9' (ubuntu)
    wgse.style.configure('TNotebook', background=def_bg)

    font = wgse.fonts.table         # Must have Fonts subsystem setup by now

    # dna frame (Program Title area above Tabs)
    dnaFrame = Frame(wgse.window, relief="solid")
    dnaFrame.grid(row=0, padx=5, pady=2)
    dnaFrame.columnconfigure(0, weight=1)
    dnaFrame.rowconfigure(0, weight=1)

    # Logo
    wgse.dnaImage = ImageTk.PhotoImage(Image.open(wgse.image_oFP))      # Needs stability outside this call
    dnaLabel = Label(dnaFrame, compound="center", text="", image=wgse.dnaImage)
    dnaLabel.grid(row=0, column=0, padx=1, pady=1)  # command=lambda: webbrowser.open_new(wgse.manual_url)

    wgse.deflab_fg = dnaLabel.cget("foreground")     # Saved default label text color

    # Title (program name, version, manual / exit buttons, current set file name)
    headlineFrame = Frame(dnaFrame)
    headlineFrame.grid(row=0, column=1, padx=1)

    Label(headlineFrame, text="WGS Extract", font=font['32']).grid(column=0, row=0, columnspan=3)
    version_fg = 'red' if is_wgse_outdated() else wgse.deflab_fg

    documentationButton = Button(headlineFrame, text=wgse.lang.i18n['WGSExtractManual'],
                                 font=font['13'], command=lambda: webbrowser.open_new(wgse.manual_url))
    documentationButton.grid(column=0, row=1)

    versionLabel = Label(headlineFrame, text=F' {wgse.__version__} ', font=font['12b'], fg=version_fg)
    versionLabel.grid(column=1, row=1)

    exitButton = Button(headlineFrame, text=wgse.lang.i18n['ExitProgram'], font=font['13'], command=button_exit)
    exitButton.grid(column=2, row=1)

    wgse.defbut_bg = exitButton.cget("background")       # Saved default button background so can revert colored buttons

    titleFileLabel = Label(headlineFrame, text="", font=font['12'])
    titleFileLabel.grid(column=0, row=2, columnspan=3, pady=3)

    titleCommandLabel = Label(headlineFrame, text="", font=font['12'])
    titleCommandLabel.grid(column=0, row=3, columnspan=3, pady=3)

    # Setup Tabs area (notebook sub-frames) for Main Window (TTK feature)
    tabParent = Notebook(wgse.window)
    tabParent.grid(row=1, padx=10, pady=5)
    tabParent.columnconfigure(0, weight=1)
    tabParent.rowconfigure(0, weight=1)
    tabParent.enable_traversal()

    # Setup each Tab as a Frame in the Parent Notebook (TTK Feature)
    tabSettings = Frame(tabParent)
    tabExtract = Frame(tabParent)
    tabAnalyze = Frame(tabParent)
    tabDEBUG = Frame(tabParent)
    tabParent.add(tabSettings, text=wgse.lang.i18n['TabSettings'])
    tabParent.add(tabExtract, text=wgse.lang.i18n['TabExtractData'])
    tabParent.add(tabAnalyze, text=wgse.lang.i18n['TabAnalyze'])
    tabParent.add(tabDEBUG, text="DEBUG")
    if not wgse.DEBUG_MODE:
        tabParent.hide(3)

    ##################################################################################################################
    # Settings Tab

    # Settings frame        (these buttons always available; so no state set)
    fileSelectFrame = LabelFrame(tabSettings, text=wgse.lang.i18n['FrameSettings'], font=font['14'])
    fileSelectFrame.grid(row=1, padx=10, pady=5, sticky=NSEW)
    fileSelectFrame.columnconfigure(1, weight=1)
    fileSelectFrame.rowconfigure(3, weight=1)
    crow = 0

    # Output Directory setting
    outdirLabel = Label(fileSelectFrame, text=wgse.lang.i18n['OutputDirectory'], font=font['14'])
    outdirLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.outdir.FB if wgse.outdir and wgse.outdir.FP and wgse.outdir.user_set else wgse.lang.i18n['Default']
    outdirButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_select_outdir)
    outdirButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W)
    outdirUnselectButton = Button(fileSelectFrame, text=wgse.lang.i18n['Unselect'], font=font['14'],
                                  command=button_unselect_outdir)
    outdirUnselectButton.grid(column=3, row=crow, padx=5, pady=2); crow += 1
    if not (wgse.outdir and wgse.outdir.oFP):
        outdirUnselectButton.grid_remove()         # Will not display if no directory selected

    # Reference Library Directory setting (default: reference_library in installation directory)
    reflibDirectoryLabel = Label(fileSelectFrame, text=wgse.lang.i18n['ReferenceLibrary'], font=font['14'])
    reflibDirectoryLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.reflib.FB if wgse.reflib.set else wgse.lang.i18n['Default']
    reflibDirectoryButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_set_referencelibrary)
    reflibDirectoryButton.grid(column=1, columnspan=3, row=crow, padx=5, pady=2, sticky=W); crow += 1

    # Temporary Files Directory setting (default: temp directory in installation directory)
    tempDirectoryLabel = Label(fileSelectFrame, text=wgse.lang.i18n['TempDirectory'], font=font['14'])
    tempDirectoryLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.tempf.FB if wgse.tempf and wgse.tempf.set else wgse.lang.i18n['Default']
    tempDirectoryButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_set_tempdir)
    tempDirectoryButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W)

    prefServerButton = Button(fileSelectFrame, text=wgse.prefserver, font=font['14'],
                              command=button_toggle_preferred_server)
    prefServerButton.grid(column=3, row=crow, padx=5, pady=5, stick=E); crow += 1

    # Language setting (with language reload and debug_mode toggle)
    langSelectLabel = Label(fileSelectFrame, text=wgse.lang.i18n["LanguageSetting"], font=font['14'])
    langSelectLabel.grid(column=0, row=crow, padx=5, pady=2)
    langSelectButton = Button(fileSelectFrame, text=wgse.lang.language, font=font['14'],
                              command=button_change_language)
    if not wgse.DEBUG_MODE:
        langSelectButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W)
    else:
        langSelectButton.grid(column=1, row=crow, padx=5, pady=2, sticky=W)
        lreloadButton = Button(fileSelectFrame, text=wgse.lang.i18n["LanguageReload"], font=font['14'],
                               command=button_language_reload)
        lreloadButton.grid(column=2, row=crow, padx=5, pady=2)

    selector = "buttonDebugModeOff" if wgse.DEBUG_MODE else "buttonDebugModeOn"
    debugToggleButton = Button(fileSelectFrame, text=wgse.lang.i18n[selector], font=font['14'],
                               command=button_debug_mode_toggle)
    debugToggleButton.grid(column=3, row=crow, padx=5, pady=2, sticky=E); crow += 1

    # BAM File Selection and Stats      (Buttons initially disabled until enabled by other actions)
    fileInfoFrame = LabelFrame(tabSettings, text=wgse.lang.i18n['FrameBAM'], font=font['14'])
    fileInfoFrame.grid(row=2, padx=10, pady=5, sticky=NSEW)
    fileInfoFrame.columnconfigure(1, weight=1)
    fileInfoFrame.rowconfigure(0, weight=1)
    crow = 0

    # We configure with one main column of labels on the left and one main column of data on the right.
    # But we further divide so the left column can have two columns of buttons and the right can have 3 columns of buttons.
    # So we columnspan and give weights to equally space them.
    fileInfoFrame.grid_columnconfigure(0, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(1, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(2, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(3, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(4, weight=1, uniform="mybuttons")  # Let last take slack
    fileLabelColSpan = 2    # Column span of left label column
    fileInfoColSpan = 3     # Column span of right data column

    # May have restored BAM from settings before calling here; or recreating window due to language change
    bamSelectedLabel = Label(fileInfoFrame, text=wgse.lang.i18n['FilenameOfSourceBam'], font=font['14'])
    bamSelectedLabel.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamSelectButton = Button(fileInfoFrame, text=wgse.lang.i18n['SelectBamFile'], font=font['14'],
                             command=button_select_BAM_file, state="disabled")
    bamSelectButton.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan, sticky=W); crow += 1

    infoLabelReferenceGenomeOfBam = Label(fileInfoFrame, text=wgse.lang.i18n['BAMAlignedToReferenceGenome'],
                                          font=font['14'])
    infoLabelReferenceGenomeOfBam.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamReferenceGenomeLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamReferenceGenomeLabel.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelMapAvgReadDepth = Label(fileInfoFrame, text=wgse.lang.i18n['MapAvgReadDepthNoN'] + ":", font=font['14'])
    infoLabelMapAvgReadDepth.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamMapAvgReadDepthLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamMapAvgReadDepthLabel.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelGender = Label(fileInfoFrame, text=wgse.lang.i18n['Gender'] + ":", font=font['14'])
    infoLabelGender.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamGenderLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamGenderLabel.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelChroms = Label(fileInfoFrame, text=wgse.lang.i18n['Chroms'] + ":", font=font['14'])
    infoLabelChroms.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamChromsLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamChromsLabel.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelFsize = Label(fileInfoFrame,
                           text=f'{wgse.lang.i18n["FileStats"]}:', font=font['14'])
    infoLabelFsize.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)
    bamFsizeLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamFsizeLabel.grid(column=2, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    # Label for buttons in BAM File Frame
    bamIdxstatsLabel = Label(fileInfoFrame, text=wgse.lang.i18n['BamCheckSamtoolsDescription'], font=font['14'])
    bamIdxstatsLabel.grid(column=0, row=crow, padx=5, pady=2, columnspan=fileLabelColSpan)

    # Realign Button below label for all buttons here
    bamRealignButton = Button(fileInfoFrame, text=wgse.lang.i18n['Realign'], font=font['14'],
                              command=button_realign_BAM, state="disabled")
    bamRealignButton.grid(column=0, row=crow+1, padx=5, pady=2)
    bamUnselectButton = Button(fileInfoFrame, text=wgse.lang.i18n['Unselect'], font=font['14'],
                               command=button_unselect_BAM_file, state="disabled")
    bamUnselectButton.grid(column=1, row=crow+1, padx=5, pady=2)
    bamRealignButton.grid_remove()
    bamUnselectButton.grid_remove()

    # Stats button is special. Normally we call stats when BAM selected and this just shows results. But if not
    #  sorted or indexed, or a CRAM file, then stats is not run by default. This button overrides that and runs
    #  stats anyway. But we encourage the user through the warning to hit the sort and/or index button or convert the
    #  CRAM to a BAM first instead of forcing the stats to run otherwise here.
    bamIdxstatsButton = Button(fileInfoFrame, text=wgse.lang.i18n['BamCheckSamtoolsButton'], font=font['14'],
                               command=button_stats_BAM, state="disabled")
    bamIdxstatsButton.grid(column=2, row=crow, padx=5, pady=2)

    # Header view button
    bamHeaderButton = Button(fileInfoFrame, text=wgse.lang.i18n['BamHeader'], font=font['14'],
                             command=button_header_BAM, state="disabled")
    bamHeaderButton.grid(column=2, row=crow+1, padx=5, pady=2)

    # Sort and Index buttons are special.  Disabled and set as done if BAM already that way.
    # Will be enabled and set as to do when BAM is not that way.
    bamIndexButton = Button(fileInfoFrame, text=wgse.lang.i18n['Indexed'], font=font['14'],
                            command=button_index_BAM, state="disabled")
    bamIndexButton.grid(column=3, row=crow, padx=5, pady=2)
    bamIndexButton.grid_remove()

    bamSortButton = Button(fileInfoFrame, text=wgse.lang.i18n['Sorted'], font=font['14'],
                           command=button_sort_BAM, state="disabled")
    bamSortButton.grid(column=3, row=crow+1, padx=5, pady=2)
    bamSortButton.grid_remove()

    # Just setup now as if a BAM file specified; will set appropriately before grid() to display
    bamConvertButton = Button(fileInfoFrame, text=wgse.lang.i18n['ToCRAM'], font=font['14'],
                              command=button_BAM_to_CRAM, state="disabled")
    bamConvertButton.grid(column=4, row=crow, padx=5, pady=2)
    bamConvertButton.grid_remove()

    bamWESButton = Button(fileInfoFrame, text=wgse.lang.i18n['ToWESBAM'], font=font['14'],
                              command=button_extract_WES, state="disabled")
    bamWESButton.grid(column=4, row=crow+1, padx=5, pady=2)
    bamWESButton.grid_remove()

    ##################################################################################################################
    # Extract Tab

    Label(tabExtract, text=wgse.lang.i18n['SuggestedFor'], font=font['14u']).grid(column=0, row=0, padx=5, pady=2)
    Label(tabExtract, text=wgse.lang.i18n['GenerateFile'], font=font['14u']).grid(column=1, row=0, padx=5, pady=2)

    # ------------ Microarray / Autosomal Frame -------------------------
    autosomesFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['FrameMicroarray'], font=font['16'])
    autosomesFrame.grid(row=2, columnspan=2, padx=10, pady=2, sticky=NSEW)
    autosomesFrame.columnconfigure(0, weight=1)
    autosomesFrame.rowconfigure(1, weight=1)

    autosomesLabel = Label(autosomesFrame, text=wgse.lang.i18n['AutosomesDescr'], font=font['14'])
    autosomesLabel.grid(column=0, row=0, rowspan=2, padx=5, pady=2)

    autosomalFormatsButton = Button(autosomesFrame, text=wgse.lang.i18n['SelectAutosomalFormats'], font=font['14'],
                                    command=button_microarray_window, state="disabled")
    autosomalFormatsButton.grid(column=1, row=0, padx=5, pady=2)

    autosomesTrailerLabel = Label(autosomesFrame, text=wgse.lang.i18n['AutosomesTrailer'], font=font['13'])
    autosomesTrailerLabel.grid(column=0, row=3, columnspan=2, padx=5, pady=2)

    # ------------ Mitochondial Frame -------------------------
    mitoFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['FrameMitochondrialDNA'], font=font['16'])
    mitoFrame.grid(row=3, columnspan=2, padx=10, pady=2, sticky=NSEW)
    mitoFrame.columnconfigure(0, weight=1)
    mitoFrame.rowconfigure(0, weight=1)

    mtLabel = Label(mitoFrame, text=wgse.lang.i18n['MitochondrialDescriptionNeededFor'], font=font['14'])
    mtLabel.grid(column=0, row=0, padx=5, pady=2)
    mitoFASTAButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateFastaMtdna'], font=font['14'],
                             command=button_mtdna_FASTA, state="disabled")
    mitoFASTAButton.grid(column=1, row=0, padx=5, pady=2)

    mtLabel2 = Label(mitoFrame, text=wgse.lang.i18n['RequiredForMitoverse'], font=font['14'])
    mtLabel2.grid(column=0, row=1, rowspan=2, padx=5, pady=2)

    mitoBAMButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateBAMOnlyWithMT'], font=font['14'],
                             command=button_mtdna_BAM, state="disabled")
    mitoBAMButton.grid(column=1, row=1, padx=5, pady=2)

    mitoVCFButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateVCFOnlyWithMT'], font=font['14'],
                             command=button_mtdna_VCF, state="disabled")
    mitoVCFButton.grid(column=1, row=2, padx=5, pady=2)

    # mtLabel2 = Label(mitoFrame, text=wgse.lang.i18n['MitochondrialDescriptionNeededFor2'], font=font['14'])
    # mtLabel2.grid(column=0, row=3, padx=5, pady=2)

    # ------------ Y Chromosome Frame -------------------------
    yFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['FrameYDNA'], font=font['16'])
    yFrame.grid(row=4, columnspan=2, padx=10, pady=2, sticky=NSEW)
    yFrame.columnconfigure(0, weight=1)
    yFrame.rowconfigure(1, weight=1)

    Label(yFrame, text=wgse.lang.i18n['RequiredForUploadingToYfull'], font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    yANDmtButton = Button(yFrame, text=wgse.lang.i18n['GenerateBAMwithYandMtdna'], font=font['14'], command=button_yAndMt_BAM,
                          state="disabled")
    yANDmtButton.grid(column=1, row=0, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForYDNAWarehouse'], font=font['14']).grid(column=0, row=1, padx=5, pady=2)
    yOnlyButton = Button(yFrame, text=wgse.lang.i18n['GenerateBAMOnlyWithY'], font=font['14'], command=button_yOnly_BAM,
                         state="disabled")
    yOnlyButton.grid(column=1, row=1, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForCladefinder'], font=font['14']).grid(column=0, row=2, padx=5, pady=2)
    yVCFButton = Button(yFrame, text=wgse.lang.i18n['GenerateVCFOnlyWithY'], font=font['14'], command=button_yOnly_VCF,
                        state="disabled")
    yVCFButton.grid(column=1, row=2, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForUploadingToYfull2'],
          font=font['13']).grid(columnspan=2, row=3, padx=5, pady=2)

    ##################################################################################################################
    # Analyze Tab

    haplogroupFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['FrameHaplogroups'], font=font['16'])
    haplogroupFrame.grid(row=2, padx=10, pady=5, sticky=NSEW)
    haplogroupFrame.columnconfigure(0, weight=1)
    haplogroupFrame.rowconfigure(1, weight=1)

    Label(haplogroupFrame, text=wgse.lang.i18n['DetermineHaplogroups'],
          font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    haplogroupYButton = Button(haplogroupFrame, text=wgse.lang.i18n['YDNA'], font=font['14'],
                               command=button_ydna_haplogroup, state="disabled")
    haplogroupYButton.grid(column=1, row=0, padx=5, pady=2)
    haplogroupMtButton = Button(haplogroupFrame, text=wgse.lang.i18n['MitochondrialDNA'], font=font['14'],
                                command=button_mtdna_haplogroup, state="disabled")
    haplogroupMtButton.grid(column=2, row=0, padx=5, pady=2)
    Label(haplogroupFrame, text=wgse.lang.i18n['DetermineHaplogroups2'],
          font=font['13']).grid(columnspan=3, row=1, padx=5, pady=2)

    oralMicrobiomeFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['FrameOralMicrobiome'], font=font['16'])
    oralMicrobiomeFrame.grid(row=3, padx=10, pady=5, sticky=NSEW)
    oralMicrobiomeFrame.columnconfigure(0, weight=1)
    oralMicrobiomeFrame.rowconfigure(1, weight=1)

    Label(oralMicrobiomeFrame, text=wgse.lang.i18n['ForUploadsCosmosId'],
          font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    exportUnmappedReadsButton = Button(oralMicrobiomeFrame, text=wgse.lang.i18n['ExportUnmappedReads'], font=font['14'],
                                       command=button_export_unmapped_reads, state="disabled")
    exportUnmappedReadsButton.grid(column=1, row=0, padx=5, pady=2)

    # FASTQ Files Frame
    fastqFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['FrameFASTQFiles'], font=font['16'])
    fastqFrame.grid(row=4, padx=10, pady=5, sticky=NSEW)
    fastqFrame.columnconfigure(0, weight=1)
    fastqFrame.rowconfigure(1, weight=1)

    bamUnalignLabel = Label(fastqFrame, text=wgse.lang.i18n['UnalignTitle'], font=font['14'])
    bamUnalignButton = Button(fastqFrame, text=wgse.lang.i18n['Unalign'], font=font['14'],
                              command=button_unalign_BAM, state="disabled")
    bamAlignButton = Button(fastqFrame, text=wgse.lang.i18n['Align'], font=font['14'],
                            command=button_align_BAM, state="disabled")
    fastqFastpButton = Button(fastqFrame, text=wgse.lang.i18n['Fastp'], font=font['14'],
                              command=button_fastp_FASTQ, state="disabled")
    fastqFastqcButton = Button(fastqFrame, text=wgse.lang.i18n['Fastqc'], font=font['14'],
                               command=button_fastqc_FASTQ, state="disabled")

    bamUnalignLabel.grid(row=0, column=0, padx=5, pady=2)
    bamUnalignButton.grid(row=0, column=2, padx=5, pady=2)
    bamAlignButton.grid(row=0, column=3, padx=5, pady=2)
    fastqFastpButton.grid(row=1, column=2, padx=5, pady=2)
    fastqFastqcButton.grid(row=1, column=3, padx=5, pady=2)

    # VCF Files Frame
    vcfFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['FrameVCFFiles'], font=font['16'])
    vcfFrame.grid(row=5, padx=10, pady=5, sticky=NSEW)
    vcfFrame.columnconfigure(0, weight=1)
    vcfFrame.rowconfigure(1, weight=1)

    vcfGenerateLabel = Label(vcfFrame, text=wgse.lang.i18n['GenerateVCFTitle'], font=font['14'])
    vcfGenerateLabel.grid(row=0, column=0, padx=5, pady=2)
    vcfModifyLabel = Label(vcfFrame, text=wgse.lang.i18n['ModifyVCFTitle'], font=font['14'])
    vcfModifyLabel.grid(row=1, column=0, padx=5, pady=2)

    SNPVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFSNP'], font=font['14'],
                          command=button_SNP_VCF, state="disabled")
    InDelVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFInDel'], font=font['14'],
                            command=button_InDel_VCF, state="disabled")
    CNVVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFCNV'], font=font['14'],
                          command=button_CNV_VCF, state="disabled")
    SVVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFSV'], font=font['14'],
                         command=button_SV_VCF, state="disabled")

    if wgse.DEBUG_MODE:
        SNPVCFButton.grid(column=1, row=0, padx=5, pady=2)
        InDelVCFButton.grid(column=2, row=0, padx=5, pady=2)
        CNVVCFButton.grid(column=3, row=0, padx=5, pady=2)
        SVVCFButton.grid(column=4, row=0, padx=5, pady=2)

    AnnotateVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFAnnotate'], font=font['14'],
                               command=button_annotate_VCF, state="disabled")
    AnnotateVCFButton.grid(column=1, columnspan=2, row=1, padx=5, pady=2)
    FilterVCFButton = Button(vcfFrame, text=wgse.lang.i18n['GenerateVCFFilter'], font=font['14'],
                             command=button_filter_VCF, state="disabled")
    FilterVCFButton.grid(column=3, row=1, padx=5, pady=2)
    VarQCButton = Button(vcfFrame, text=wgse.lang.i18n['VariantQC'], font=font['14'],
                         command=button_varqc_VCF, state="disabled")
    VarQCButton.grid(column=4, row=1, padx=5, pady=2)

    # Todo "Add to" check box option so it does not create a newly named VCF but instead modifies in-place

    ##################################################################################################################
    # DEBUG Tab

    # Debug Frame will only appear (get gridded) if DEBUG_MODE is on; but still need to create first here
    debugFrame = LabelFrame(tabDEBUG, text=wgse.lang.i18n['FrameDEBUG_MODE'], font=font['14'])
    debugFrame.columnconfigure(0, weight=1)
    debugFrame.rowconfigure(1, weight=1)

    # Button to subset a BAM by a percentage
    subsetBAMLabel = Label(debugFrame, text=wgse.lang.i18n["SubsetBAM"], font=font['14'])
    subsetBAMButton = Button(debugFrame, text=wgse.lang.i18n['buttonSubsetBAM'],
                             font=font['14'], command=_button_SubsetBAM)

    # Button to run microarray generation in parallel (experimental; known problems with hsxx models)
    runMicroParallelLabel = Label(debugFrame, text=wgse.lang.i18n['buttonAllSNPsParallel'], font=font['14'])
    runMicroParallelButton = Button(debugFrame, text=wgse.lang.i18n['buttonGenerateSelectedFiles'], font=font['14'],
                                    command=lambda win=wgse.window: _button_AllSNPs(win, parallel=True))

    # Windows 10 WSL2 BWA Aligner Override
    wsllabel = wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else wgse.lang.i18n['Inactive']
    wslbwaLabel = Label(debugFrame, text=wgse.lang.i18n['Win10WslBwaOverride'], font=font['14'])
    wslbwaButton = Button(debugFrame, text=wsllabel, font=font['14'], command=button_wsl_bwa_patch)

    # Buttons to set user overrides of total memory and threads available to WGSE
    maxsettingsLabel = Label(debugFrame, text=wgse.lang.i18n["MaxSettings"], font=font['14'])

    maxmemButton = Button(debugFrame, text=str(wgse.os_totmem_saved),
                          font=font['14'], command=_button_maxmem_setting)
    tlabel = f'{wgse.lang.i18n["MaxMem"]} {wgse.lang.i18n["Of"]} {wgse.os_totmem_proc//10**9},'
    maxmemLabel = Label(debugFrame, text=tlabel, font=font['14'])

    maxthreadButton = Button(debugFrame, text=str(wgse.os_threads_saved),
                             font=font['14'], command=_button_maxthread_setting)
    tlabel = f'{wgse.lang.i18n["MaxThread"]} {wgse.lang.i18n["Of"]} {wgse.os_threads_proc}'
    maxthreadLabel = Label(debugFrame, text=tlabel, font=font['14'])

    # Buttons to set user overrides for base font point size and typeface
    curfont=(wgse.fonts.face, wgse.fonts.basept+2)      # Nominally use 14 (base 12 + 2)
    fontsetLabel = Label(debugFrame, text=wgse.lang.i18n["SetBaseFont"], font=curfont)
    fontsizeButton = Button(debugFrame, text=str(wgse.fonts.basept)+' pt', font=curfont, command=_button_fontsize)
    flabel = wgse.fonts.face[:22]+'...' if len(wgse.fonts.face) > 25 else wgse.fonts.face
    fontfaceButton = Button(debugFrame, text=flabel, font=curfont, command=_button_fontface)
    wresetButton = Button(debugFrame, text=wgse.lang.i18n["WindowRedraw"], font=curfont,
                          command=mainwindow_reset)

    # Ugly; for translators only. Static language selection buttons in case they really screw up the current language
    lselLabel = Label(debugFrame, text=wgse.lang.i18n['RequestLanguage'], font=font['14'])
    lselButton = [x for x in range(len(wgse.lang.avail_langs))]
    langcol = 0
    for lang in wgse.lang.avail_langs:  # Need all languages so go back to original _lang dictionary
        # lselLabel[column] = Label(debugFrame, text=wgse.lang._lang[lang]["RequestLanguage"], font=font['14'])
        # noinspection PyTypeChecker
        lselButton[langcol] = Button(debugFrame, text=wgse.lang.avail_langs[langcol], font=font['14'],
                                    command=lambda l=lang: wgse.lang.switch_language(l))
        langcol += 1

    releaseLabel = Label(debugFrame, text=wgse.lang.i18n["ReleaseTrack"], font=font['14'])
    releaseTrackButton = Button(debugFrame, text=wgse.lang.i18n[wgse.track], font=font['14'],
                                command=_button_change_release_track)
    releaseUpdateReleaseButton = Button(debugFrame, text=wgse.lang.i18n['UpdateRelease'], font=font['14'],
                                        command=update_release)

    # debugFrame and content only becomes visible when we grid it; so do that now if debug mode is on
    if wgse.DEBUG_MODE:
        debugFrame.grid(row=6, padx=10, pady=5, sticky=NSEW)

        crow=0
        subsetBAMLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        subsetBAMButton.grid(row=crow, column=3, padx=5, pady=2, sticky=W)

        crow += 1
        runMicroParallelLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        runMicroParallelButton.grid(row=crow, column=3, padx=5, pady=2, sticky=W)

        crow += 1
        wslbwaLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        wslbwaButton.grid(row=crow, column=3, padx=5, pady=2, sticky=W)

        crow += 1
        maxsettingsLabel.grid(row=crow, column=0, padx=5, pady=2)
        maxmemButton.grid(row=crow, column=1, padx=5, pady=2, sticky=E)
        maxmemLabel.grid(row=crow, column=2, padx=5, pady=2, sticky=W)
        maxthreadButton.grid(row=crow, column=3, padx=5, pady=2, sticky=E)
        maxthreadLabel.grid(row=crow, column=4, padx=5,   pady=2, sticky=W)

        crow += 1
        fontsetLabel.grid(row=crow, column=0, padx=5, pady=2)
        fontsizeButton.grid(row=crow, column=1, padx=5, pady=2, sticky=E)
        fontfaceButton.grid(row=crow, column=2, columnspan=2, padx=5, pady=2, sticky=W)
        wresetButton.grid(row=crow, column=4,  padx=5, pady=2)

        # Create rows of 4 buttons each for various languages; large lanugage names span columns
        crow += 1  ;  ccol = 0
        lselLabel.grid(row=crow, column=0, padx=5, pady=2)
        maxchar = 20  ;  maxcol = 4     # Each column is max 20 characters; max of 5 columns
        for numlangs in range(len(wgse.lang.avail_langs)):      # For each language selectionbutton
            # Figure out how many columns a button needs
            # noinspection PyUnresolvedReferences
            span = min((len(lselButton[numlangs].cget("text")) // maxchar) + 1, maxcol)
            if ccol + span > maxcol:     # Not enough room in row left for this button; start a new row
                crow += 1
                ccol = 0
            # noinspection PyUnresolvedReferences
            lselButton[numlangs].grid(row=crow, column=ccol+1, columnspan=span, padx=5, pady=2)
            ccol = (ccol + span) % maxcol
            crow += 1 if ccol == 0 else 0

        crow += 1
        releaseLabel.grid(row=crow, column=0, padx=5, pady=2)
        releaseTrackButton.grid(row=crow, column=1, padx=5, pady=2)
        releaseUpdateReleaseButton.grid(row=crow, column=2, padx=5, pady=2)

    # Todo add theme setting and button to change it  (in development; may change GUI from tkinter before completing)
    ''' Developing code to allow a theme change setting
    style = ttk.Style()
    cur_theme = style.theme_use()
    themes = StringVar(value=style.named_themes())

    Label(debugFrame, text=f"Change Theme (current highlighted):", font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    select_theme = Listbox(debugFrame, listvariable=themes, selectmode="browse", height=1)
    select_theme.grid(column=1, row=0, padx=5, pady=2)
    select_theme.selection_set(themes.index[cur_theme])
    select_theme.bind('<<ListboxSelect>>', change_theme)
    '''

    # Colorize long run buttons; pink for an one to three hours; deep red for 4+ hours (Align)
    # Todo Do not color if results file available or specific direction that is less (e.g. CRAM to BAM)
    bamRealignButton.config(bg=red)
    bamAlignButton.config(bg=red)
    autosomalFormatsButton.config(bg=pink)
    fastqFastqcButton.config(bg=pink)
    bamConvertButton.config(bg=pink)
    SNPVCFButton.config(bg=pink)
    InDelVCFButton.config(bg=pink)
    CNVVCFButton.config(bg=pink)
    SVVCFButton.config(bg=pink)

    # If default value has been changed (overriden) by user (setting), then color button to indicate so
    override = 2 <= wgse.os_totmem_saved < wgse.os_totmem_proc // 10**9
    maxmemButton.configure(bg=cyan if override else wgse.defbut_bg)
    override = 1 <= wgse.os_threads_saved < wgse.os_threads_proc
    maxthreadButton.configure(bg=cyan if override else wgse.defbut_bg)
    override = wgse.fonts.default_basept != wgse.fonts.basept
    fontsizeButton.configure(bg=cyan if override else wgse.defbut_bg)
    override = wgse.fonts.default_face != wgse.fonts.face
    fontfaceButton.configure(bg=cyan if override else wgse.defbut_bg)

    # Listing of main window action buttons that may need the state changed based on content of the loaded BAM file
    # Does not include Settings buttons nor BAM file button itself which are handled separately.
    all_action_buttons = [
        # Settings Tab (inside BAM File setting only)
        # bamSelectButton,       # Not a typical action button. Always available unless outputdir not set.
        bamIdxstatsButton,
        bamHeaderButton,
        bamIndexButton,
        bamSortButton,
        bamConvertButton,
        bamWESButton,
        bamRealignButton,
        bamUnselectButton,
        # Extract Data Tab
        autosomalFormatsButton,
        mitoFASTAButton,
        mitoBAMButton,
        mitoVCFButton,
        yANDmtButton,
        yOnlyButton,
        yVCFButton,
        # Analyze Tab
        haplogroupYButton,
        haplogroupMtButton,
        exportUnmappedReadsButton,
        bamUnalignButton,
        # bamAlignButton,       # Not a typical action button. Always available unless outputdir not set.
        # fastqFastpButton,
        # fastqFastqcButton,
        SNPVCFButton,
        InDelVCFButton,
        CNVVCFButton,
        SVVCFButton,
        # AnnotateVCFButton,     # Not a typical action button. Always available unless outputdir not set.
        # FilterVCFButton,       # Not a typical action button. Always available unless outputdir not set.
        # VarQCButton,           # Not a typical action button. Always available unless outputdir not set.
        # DEBUG Tab
        # lselButton[numlangs]   # Language buttons are dynamic and dependent on content of languages.xlsx file
        # wslbwaButton,
        runMicroParallelButton,
        subsetBAMButton
        # maxmemButton,
        # maxthreadButton,
        # fontsizeButton,
        # fontfaceButton,
        # wresetButton
    ]

    # We may have called mainwindow_setup() after already loading saved settings (including a BAMFile) or
    #   due to a language change midstream. So need to modify default setup for changed values.
    set_BAM_window_settings()      # If BAM not set from stored settings; simply makes sure settings are cleared
    mainwindow_resume()            # Includes call to update_action_buttons()
    return wgse.window
