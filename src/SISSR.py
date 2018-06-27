'''
@attention: using SISSR peak finder in PFMeraserver for finding peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 18 Jan 2011
@version: 1.0
@todo: produce wiggle output format as well

'''

import os
from src.FunctionsUtility import print_message

def SISSR_call(label, data_path, pfms_path, wiggle, control_file, sissr_options):

    language = "perl"
    space = " "

    print_message('SISSr is Starting')

    if control_file != "":
        control_file = "-b " + control_file
    if wiggle == "-wig":
        print "-wig is not supported for SISSR"
        '''SISSR_call_command = language + space + SISSR_path + space +  "-i" + space + data_path + space +"-o" + space + "sissr_" + label + ".bed" + space + control_file + space + sissr_options        
            os.system(SISSR_call_command)'''
    else:
        SISSR_call_command = language + space + pfms_path + space +  "-i" + space + data_path + space +"-o" + space + label + "_sissr" + ".bed" + space + control_file + space + sissr_options
        
        os.system(SISSR_call_command)
        print "check the output file"
        
        SISSR_output_file = label + "_sissr" + ".bed"
    
    print_message('SISSr is Exiting')
    return SISSR_output_file
