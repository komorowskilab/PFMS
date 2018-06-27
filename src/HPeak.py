'''
Created on May 29, 2010
@author: kruczyk

@attention: using HPeak software in PFMeraserver for identifying peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 19 Dec 2010
@version: 1.0
@note: HPeak is giving divison by zero error when control data was used.

'''

import os
import shutil
from src.FunctionsUtility import print_message
from src import FunctionsUtility

def HPeak_data_file_parser(path):
    HPeak_parsed_data_file = "./HPeak_parsed_data_file.bed"
    BED_data_file_id = open(path, "r")
    HPeak_parsed_data_file_id = open(HPeak_parsed_data_file, "w")
    line = BED_data_file_id.readline()
    output_string = " "
    while line:
        list = line.split('\t')
        output_string = list[0] + '\t' + list[1] + '\t' + list[2] + '\t' + list[5]
        HPeak_parsed_data_file_id.write(output_string)
        line = BED_data_file_id.readline()
    BED_data_file_id.close()
    HPeak_parsed_data_file_id.close()
    return HPeak_parsed_data_file
    
def HPeak_call(label, data_path, pfms_path, wiggle, control_path, hpeak_options):
    control_path = ""
    print_message('HPeak is Starting')
    HPeak_format = " -format BED"
    HPeak_option_t = " -t"
    HPeak_data_path =  "tmp.dat"
    HPeak_option_n = " -n"
    perl = "perl"
    space = " "
    if wiggle == "-wig":
        HPeak_option_wig = " -wig"
    else:
        HPeak_option_wig = ""
    
    #parse the data file      
    parsed = HPeak_data_file_parser(data_path)
    full_data_path = "./" + parsed
    full_data_path = os.path.expanduser(full_data_path)
    HPeak_create_tmp_file_command = "echo " + full_data_path + ">" + HPeak_data_path
    os.system(HPeak_create_tmp_file_command)
    #parse the control file
    HPeak_control_option = ""
    parsed_control = ""
    if control_path != "":
        parsed_control = HPeak_data_file_parser(control_path)
        HPeak_control_path = "tmp_control.dat"
        full_control_path = "./" + parsed_control
        full_control_path = os.path.expanduser(full_control_path)
        HPeak_create_tmp_control_file_command = "echo " + full_control_path + ">" + HPeak_control_path
        HPeak_control_option = "-c " + HPeak_control_path
        os.system(HPeak_create_tmp_control_file_command)
        
    HPeak_call_command = perl + space + pfms_path + HPeak_format + HPeak_option_t + space + HPeak_data_path + HPeak_option_n + space + label + "_HPeak " + space + HPeak_option_wig + space + hpeak_options
   
#    HPeak_call_command = perl + space + pfms_path + HPeak_format + HPeak_option_t + space + HPeak_data_path + HPeak_option_n + space + label + "_HPeak " + space + HPeak_control_option + space + HPeak_option_wig + space + hpeak_options
    os.system(HPeak_call_command)
    os.remove(HPeak_data_path)
    os.remove(label + "_HPeak.sum")
    os.remove(label + "_HPeak.log")
    os.remove(parsed)
    if os.path.exists(parsed_control):
        os.remove(parsed_control)
    if os.path.exists("tmp_control.dat"):
        os.remove("tmp_control.dat")
        
    if wiggle != "-wig":
        HPeak_output_file = label + "_hpeak.bed"
        try: 
            shutil.copyfile(label + "_HPeak.allregions.txt", HPeak_output_file)
        except (IOError):
            outfile = open(HPeak_output_file, 'w')
            outfile.close()
            
        if os.access(label + "_HPeak.allregions.txt", os.F_OK):
            os.remove(label + "_HPeak.allregions.txt")
        #Append chr to all the lines and chr23 to chrX, chr24 to chrY
        functionsUtility = FunctionsUtility.FunctionsUtility()
        if os.access(HPeak_output_file, os.F_OK):
            functionsUtility.append_chr(HPeak_output_file)


    else:
        HPeak_output_file = "./" + label + "_hpeak.wig"
        shutil.copyfile(label + "_HPeak.wig", HPeak_output_file)
        HPeak_regions_file = "./" + label + "_HPeak.allregions.txt"
        if os.access(label + "_HPeak.allregions.txt", os.F_OK):
            os.remove(HPeak_regions_file)
        if os.access(label + "_HPeak.wig", os.F_OK):
            os.remove(label + "_HPeak.wig")

    print_message('HPeak is Exiting')


    return HPeak_output_file
