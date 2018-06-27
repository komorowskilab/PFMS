'''
Created on May 6, 2010
@author: kruczyk

@attention: using MACS peak finder in PFMeraserver for identifying peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 1 Dec 2010
@version: 1.0

'''

import os
import shutil
from src.FunctionsUtility import print_message

def MACS_call(label, data_path, pfms_path, wiggle, control_file, macs_options):

    space = " "
    print_message('MACS is Starting')

    if control_file != "":
        control_file = "-c " + control_file

    if wiggle == "-wig":
        MACS_call_command = pfms_path +" --format BED -t " + data_path + space + control_file + space + "--name=" + label + "_macs --wig" + space + macs_options
        MACS_wig_path = "./" + label + "_macs_MACS_wiggle/treat/"
        MACS_wig_files = MACS_wig_path + label + "_macs_treat_afterfiting*"
        MACS_gunzip_command = "gunzip " + MACS_wig_files
        MACS_copy_command = "cp " + MACS_wig_files + " ./"
        os.system(MACS_call_command)
        os.system(MACS_gunzip_command)
        os.system(MACS_copy_command)
        os.system("rm " + MACS_wig_path + "*")
        os.system("rm -r ./" + label + "_macs_MACS_wiggle/")
        MACS_temp_file = "./MACS_temp_file.dat"

        MACS_echo_output_file_command = "echo " + label + "_macs_treat_afterfiting*" + ">" + MACS_temp_file
        os.system(MACS_echo_output_file_command)
        MACS_temp_file_id = open(MACS_temp_file, "r")
        MACS_output = MACS_temp_file_id.readline()
        if MACS_output[len(MACS_output) - 1:len(MACS_output)] == "\n":
            MACS_output= MACS_output[0 : len(MACS_output) - 1]
        MACS_temp_file_id.close()
        os.remove(MACS_temp_file)
        remove_macs_bed_file_command = "rm ./" + label + "_macs_peaks.bed"
        os.system(remove_macs_bed_file_command)

        MACS_output_file = label+"_macs.wig"
        os.system("cp "+ MACS_output +" "+ MACS_output_file)
        os.system("rm "+MACS_output)
    else:
        MACS_call_command = pfms_path +" --format BED -t " + data_path + space + control_file + space + "--name=" + label + "_macs" + space + macs_options
        os.system(MACS_call_command)
        MACS_output = label + "_macs_peaks.bed"
        MACS_output_file = label + "_macs.bed"
        if os.access(MACS_output, os.F_OK):
            shutil.copy(MACS_output, MACS_output_file)
            os.remove(label+"_macs_peaks.bed")

    #Removing the extra files
    if os.access(label+"_macs_negative_peaks.xls", os.F_OK):
        os.remove(label+"_macs_negative_peaks.xls")
    if os.access(label + "_macs_model.r", os.F_OK):
        os.remove(label + "_macs_model.r")
    if os.access(label + "_macs_peaks.xls", os.F_OK):
        os.remove(label + "_macs_peaks.xls")

    print_message('MACS is Exiting')
    return MACS_output_file


