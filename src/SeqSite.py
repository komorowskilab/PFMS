'''

@attention: using SeqSite peak finder in PFMeraserver for finding peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 26 Jan 2011
@version: 1.0
'''

import os
from src.FunctionsUtility import print_message

def SeqSite_call(label, data_path, pfms_path, wiggle, control_file, seqsite_options):

    language = ""
    space = " "
    
    print_message('SeqSite is Starting')

    if control_file != "":
        control_file= "-c " + control_file
    SeqSite_call_command = language + pfms_path +"SeqSite" + space +  data_path + space + label + "_seqsite" + ".bar" + space + label + "_seqsite" + ".bed" + space + control_file + space + seqsite_options
    os.system(SeqSite_call_command) #returns bed and bar files
    
    SeqSite_output_file = label + "_seqsite" + ".bed"
    
    if wiggle == "-wig":
        #conver the bar file to wig and return it
        print "WARNING: SeqSite is not working sufficiently with WIG format since size of the PFMS output would be very big when SeqSite is included."
        bar_2_wig = "perl" + space + pfms_path +"bar2wig.pl" + space + label + "_seqsite" + ".bar" + space + label + "_seqsite" + ".wig"
        os.system(bar_2_wig) #returns a wig file
        os.system("rm " + label + "_seqsite" + ".bed")
        SeqSite_output_file = label + "_seqsite" + ".wig"
        print "the bar file converted to wig"
    
    os.system("rm " + label + "_seqsite" + ".bar")
    print_message('SeqSite is Exiting')
    return SeqSite_output_file
