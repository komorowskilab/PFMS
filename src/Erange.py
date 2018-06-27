'''
Created on May 6, 2010
@author: kruczyk

@attention: using Erange peak finder in PFMeraserver for finding peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 7 Dec 2010
@version: 1.0
@todo: produce wiggle format as well, by converting the bed result file to wig format.
        or adding -wig option if supported in comming versions of Erange
'''

import os
from src.FunctionsUtility import print_message

def Erange_call(label, data_path, pfms_path, wiggle, control_file, erange_options):
    #Convert the BED input files to rds
    print_message("Erange is Starting")
    control_option = ""
    python = "python2.6 "
    Erange_rds_file = label + ".rds"
    space = " "
    Erange_parsing_command = python + pfms_path + "makerdsfrombed.py" + space + label + space + data_path + space + Erange_rds_file
    os.system(Erange_parsing_command)
    if control_file != "":
        Erange_control_rds_file =  label + "_control.rds"
        Erange_parsing_control_command = python + pfms_path +"makerdsfrombed.py"+ space + label+"_control" + space + control_file + space + Erange_control_rds_file
        os.system(Erange_parsing_control_command)
        control_option = "-control " + Erange_control_rds_file
    
    Erange_regions_file = label + "_Erange_regions.txt"
    Erange_regions_call_command = python + pfms_path +"findall.py"+space+ label + space + Erange_rds_file + space + Erange_regions_file + space + control_option +" -listPeaks" + space + erange_options
    os.system(Erange_regions_call_command)

    Erange_output_file = label + "_erange.bed"
    Erange_output_BED_command = python + pfms_path +"regiontobed.py"+ space+ label + space + Erange_regions_file + space + Erange_output_file
    os.system(Erange_output_BED_command)

    if wiggle == "-wig":
        #convert the bed to rds, then make wiggle out of the rds
        print "-wig is not supported for Erange"
        '''print "the original file has been converted to wig rather than identifying the peaks"
        wig_conversion_path = pfms_path+"makewiggle.py "    #makes wiggle file from rds file (can't recognize the result of makerdsfrombed.py as a correct rds file)
        Erange_wig_output_file = label + "_Erange.wig"
        
        old_result = python + wig_conversion_path + label + space + Erange_rds_file + space + Erange_wig_output_file
        os.system(old_result)
        #my version
        convert_bed2rds = python + pfms_path +"makerdsfrombed.py"+ space + "rdsfrombed" + space + Erange_output_file + space + "rds_wigle.rds"
        convert_rds2wig = python + wig_conversion_path + label + space + "rds_wigle.rds" + space + Erange_wig_output_file + space + erange_options
        os.system(convert_bed2rds)
        os.system(convert_rds2wig)
        print process_name, ' (Erange) is Exiting'
        return Erange_wig_output_file'''
    
    if os.path.exists(label + "_control.rds.log"):
        os.remove(label + "_control.rds.log")
        os.remove(label + "_control.rds")

    os.remove("findall.log")
    os.remove(label + "_Erange_regions.txt")
    os.remove(label+".rds.log")
    os.remove(label + ".rds")

    print_message("Erange is Exiting")
    return Erange_output_file
