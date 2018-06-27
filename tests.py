#! /usr/bin/python

# Test module

__author__="husen"
__date__ ="$Mar 24, 2011 3:34:09 PM$"

import os
import sys

def test_normal_install():
    #os.system("python setup.py install normal")
    cisgenome = False
    seqsite = False
    macs = False
    #check cisgenome installed
    if os.path.exists("./PeakFinders/cisGenome-2.0/bin/hts_peakdetectorv2"):
        cisgenome = True
    #check seqsite installed
    if os.path.exists("./PeakFinders/SeqSite1.0/SeqSite"):
        seqsite = True
    #check if macs/bin is executsble
    if os.access("./PeakFinders/MACS-1.3.7.1/bin/macs", os.X_OK):
        macs = True
    if cisgenome and seqsite and macs:
        print "The peakfinders installed under the current directory"
        return True
    if not cisgenome:
        print "CisGenome installation failed"
    if not seqsite:
        print "SeqSite installation failed"
    if not macs:
        print "MACS installation failed- you don't have sufficient permission to make it executable"


def test_install():
    #os.system("sudo python setup.py install")
    cisgenome = False
    seqsite = False
    macs = False
    peakfinders_copied = False
    #check cisgenome installed
    if os.path.exists(sys.prefix+"/PeakFinders"):
        peakfinders_copied = True
        if os.path.exists(sys.prefix+"/PeakFinders/cisGenome-2.0/bin/hts_peakdetectorv2"):
            cisgenome = True
        #check seqsite installed
        if os.path.exists(sys.prefix+"/PeakFinders/SeqSite1.0/SeqSite"):
            seqsite = True
        #check if macs/bin is executsble
        if os.access(sys.prefix+"/PeakFinders/MACS-1.3.7.1/bin/macs", os.X_OK):
            macs = True
        if cisgenome and seqsite and macs and peakfinders_copied:
            print "The peakfinders installed under "+ sys.prefix +"/PeakFinders"
            return True

    if peakfinders_copied:
        if not cisgenome:
            print "CisGenome installation failed"
            return False
        if not seqsite:
            print "SeqSite installation failed"
            return False
        if not macs:
            print "MACS installation failed- you don't have sufficient permission to make it executable"
            return False
        return True
    else:
        print "Peakfinders are not installed, please make sure you are authorized to access: "+sys.prefix
        return False

def run_test_case(message,input, output,command="PFMetaserver",control="", options="", peakfinders="-macs -cisgenome", wig=" -wig" ):
    print message +" ....Is Starting"
    os.system(command + " -i " + input +" -o "+ output +" "+control+" "+options + ""+ " "+peakfinders + " "+ wig)
    if wig =="-wig":format = "wig"
    else: format="bed"
    if os.path.exists("./Results_"+output+"/"+output+"_Results."+format):
        return True
    else:
        return False

def print_status(test_number, message):
    if test_number:
        print "Success -- "+ message
    else:
        print "Failed -- "+message

def test_cases(command, macs = "-macs ",
                cisgenome = "-cisgenome " ):
    if "no-macs" in sys.argv: macs =""
    if "no-hpeak" in sys.argv: hpeak =""
    if "no-findpeaks" in sys.argv: findpeaks =""
    if "no-seqsite" in sys.argv: seqsite=""
    if "no-sissr" in sys.argv: sissr=""
    if "no-cisgenome" in sys.argv: cisgenome=""
    if "no-erange" in sys.argv: erange=""
#    test1=run_test_case("Test1: BED format, Peak-finders: "+macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange,
#    sys.argv[1],sys.argv[2]+"-1", command, peakfinders=macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange)

    test2=run_test_case("Test2: WIG format, Peak-finders: "+macs+cisgenome,
        "Samples_without_H3K4me3.bed","wig_test11", command, peakfinders=macs+cisgenome, wig="-wig")

#    test3=run_test_case("Test3: BED format, Control data, Peak-finders: "+macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange,
#        sys.argv[1],sys.argv[2]+"-3",command, control="-control "+sys.argv[3],options="-min_rank 2", peakfinders=macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange)

#    test4=run_test_case("Test4: WIG format, Control data, Peak-finders: "+macs+cisgenome+findpeaks+hpeak,
#        sys.argv[1],sys.argv[2]+"-4", command, control="-control "+sys.argv[3], peakfinders=macs+cisgenome+findpeaks+hpeak, wig="-wig")

    #print msgs
    print "Test Summaries:"
 #   print_status(test1, "Test1: BED format, Peak-finders: "+macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange)
    print_status(test2, "Test2: WIG format, Peak-finders: "+macs+cisgenome)
  #  print_status(test3, "Test3: BED format, Control data, Peak-finders: "+macs+cisgenome+findpeaks+hpeak+seqsite+sissr+erange)
   # print_status(test4, "Test4: WIG format, Control data, Peak-finders: "+macs+cisgenome+findpeaks+hpeak)

if __name__ == "__main__":
    #usage: python tests.py input_file_path.bed output_label control_file_path.bed [-super]
    #       python tests.py [install-normal | install-super]
    if "install-normal" in sys.argv:
        test_normal_install()
    elif "install-super" in sys.argv:
        test_install()
    elif "install" not in sys.argv:
        print "Testing..."
        if "-super" in sys.argv:
            test_cases("PFMetaserver")
        else:
            test_cases("python PFMetaserver.py")

    if not float(sys.version[:3])>=2.6:
        print "PFMS doesn't work in parallel mode, please upgrade your python version to 2.6 or higher"
