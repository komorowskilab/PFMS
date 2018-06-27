'''
@attention: using CISGenome peak finder in PFMeraserver for identifying peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 14 Jan 2011
@version: 1.0

'''

import os 
from FunctionsUtility import print_message
from FunctionsUtility import get_cutoff_value
from FunctionsUtility import conver_float_to_int
from FunctionsUtility import sort_bed_file

def CISGenome_call(label, data_path, pfms_path, wiggle, control_file, cisgenome_options):

    language = ""
    space = " "
    print_message("CisGenome is Starting")

    bed_output_file = ""
    wig_output_file = ""
    if control_file != "": #Analyzing two-sample ChIP-seq data (CisGenome v2)
        #Convert input data to aln, since ALN file is the standard input for CisGenome ChIP-seq analysis 
        convertBedToAln = language + pfms_path + "bin/file_bed2aln" + space +  "-i" + space + data_path + space +"-o" + space + "CISGenome_" + label + ".aln"
        convertControlBedToAln = language + pfms_path + "bin/file_bed2aln" + space +  "-i" + space + control_file + space +"-o" + space + "CISGenome_control_" + label + ".aln"
        
        #wrtie the file paths to a tab dilimetted file, The first column gives the full file path. The second column indicates the DNA sample type (1 = IP, 0 = control).
        input_file_paths = open("cisgenome_data_path.txt", "w")
        input_file_paths.write("CISGenome_" + label + ".aln" + "\t1" + "\n")
        input_file_paths.write("CISGenome_control_" + label + ".aln" + "\t0" + "\n")
        input_file_paths.close()
        
        CisGenome_call_command =  language + pfms_path + "bin/seqpeak" + space + "-i" + space + "cisgenome_data_path.txt" + space + "-d . -o" + space + "CISGenome_" + label + space + cisgenome_options
        os.system(convertBedToAln)
        os.system(convertControlBedToAln)
        os.system(CisGenome_call_command)
        
        bed_output_file = "CISGenome_" + label + "_peak.cod"
        if os.path.exists("CISGenome_" + label + "_t.bar"):
            wig_output_file =  "CISGenome_" + label + "_t.bar"
        
        os.system("rm ./CISGenome_control_" + label + ".aln")
        os.system("rm ./cisgenome_data_path.txt")
        
    else:   # Analyzing one-sample ChIP-seq data
        #Convert input data to aln, since ALN file is the standard input for CisGenome ChIP-seq analysis 
        convertBedToAln = language + pfms_path + "bin/file_bed2aln" + space +  "-i" + space + data_path + space +"-o" + space + "CISGenome_" + label + ".aln"
        #Sort converted (.aln) data file
        sortALnDataFile = language + pfms_path + "bin/tablesorter_str" + space  + "CISGenome_" + label + ".aln" + space + "CISGenome_" + label + "_sorted" + ".aln" + ".sort"
        #convert the sorted data file to bar, this step will generate 3 bar files (_F.bar, _R.bar and .bar) but only need the .bar file in later steps
        
        convertAlnToBar = language + pfms_path + "bin/hts_aln2barv2" + space +  "-i" + space + "CISGenome_" + label + "_sorted" + ".aln" + ".sort" + space +"-o" + space + "CISGenome_" + label + ".bar"
        data_file = "CISGenome_" + label + ".bar"
        #compute FDR and choose cutoff
        fdrOptions = "-g" + space + pfms_path+"datatable/chrlist.txt" + space + "-l" + space + pfms_path+"datatable/chrlen.txt" + space + "-w" + space + "100"
        
        computeFDR = language + pfms_path + "bin/hts_windowsummaryv2" + space +  "-i" + space + data_file + space + fdrOptions + space + "-o" + space + "CISGenome_" + label + "summarys" + ".txt"
        #get fdr value from summary file
        os.system(convertBedToAln)
        os.system(sortALnDataFile)
        os.system(convertAlnToBar)
        os.system(computeFDR)
        cutoff = get_cutoff_value("CISGenome_" + label + "summarys" + ".txt")
        print 'cutoff = ', cutoff 
        #first round of peak detection, output is saved in cis_genome.cod, file (set -c to the read count level that satisfies FDR requirement. FDR for each read count level is given in the last column of summarys.txt obtained from the last step).
        peakDetectionOptions = "-c %d" % cutoff + space + cisgenome_options #" -w 100 -s 25 -br 1 -brl 30 -ssf 1" #"-w" + space + "100" + space + "-s" + space + "25" + space + "-c" + space + "10" + space + "-br" + space + "1" + space + "-brl" + space + "30" + space + "-ssf" + space + "1"
        CisGenome_call_command =  language + pfms_path + "bin/hts_peakdetectorv2" + space + "-i" + space + data_file + space + "-d ./" + space + "-o" + space +  "CISGenome_" + label  + space + peakDetectionOptions
        bed_output_file = "CISGenome_" + label + ".cod"
        wig_output_file = "CISGenome_" + label + ".bar"
        os.system(CisGenome_call_command)
    
    #check if -wig was not set then convert the .cod file to .bed    
    if wiggle == "-wig":
        if os.path.exists(wig_output_file):
            getWigOutputFile = language + pfms_path + "bin/affy_bar2wig" + space +  "-i" + space + wig_output_file + space + "-o" + space +"CISGenome_" + label +".wig"
            os.system(getWigOutputFile)
            # change last col from float to int in "CISGenome_" + label +".wig" and put write the contents into "CISGenome_wig_" + label +".wig"
            CISGenome_output_file = label + "_cisgenome" + ".wig"
            conver_float_to_int("CISGenome_" + label +".wig", CISGenome_output_file)
            print "check the output file"
            #remove extra files
            remove_extra_files = "CISGenome_" + label + "*"
            os.system("rm ./"+remove_extra_files)
            print '***************CisGenome Finished***************'
            return CISGenome_output_file
        else:
            print "CisGenome failed to create the CISGenome_" + label + "_t.bar file,\n *SeqPeak can't be used to produce wig files for selected samples."    
    else:
        getBedOutputFile = language + pfms_path + "bin/file_cod2bed" + space +  "-i" + space + bed_output_file + space + "-o" + space + label + "_cisgenome" + ".bed"
        os.system(getBedOutputFile)
        sort_bed_file(label + "_cisgenome" + ".bed", label + "_cisgenome" + ".bed")
        print "check the output file"
        #remove extra files
        remove_extra_files = "CISGenome_" + label + "*"
        os.system("rm ./"+remove_extra_files)
        
        CISGenome_output_file = label + "_cisgenome" + ".bed"
        print_message('CisGenome is Finished')
        return CISGenome_output_file
        
