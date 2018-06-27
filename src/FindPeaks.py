'''
@attention: using FindPeaks peak finder in PFMeraserver for finding peaks
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 17 Jan 2011
@version: 1.0
@note: using -bed option produces bedgraph file not bed.

'''

import os
import shutil
from src.FunctionsUtility import print_message

def FindPeaks_call(label, data_path, pfms_path, wiggle, control_file, findpeaks_options, java_max_memory_usage):

    print_message('FindPeaks is Starting')
    language = "java -jar" + java_max_memory_usage
    space = " "

    #data pre_proccessing
    control_file_fp = "control_data_fp.bed"
    control_option = ""
    if control_file != "":
        print "control data is used"
        shutil.copyfile(control_file, control_file_fp)
        sort_control_file = language + space + pfms_path +"SortFiles.jar bed ./" + space + control_file_fp
        gunzip_sorted_control_file = "gunzip -d " + control_file_fp + "*.gz"
        os.system(sort_control_file)
        os.remove(control_file_fp)
        os.system(gunzip_sorted_control_file)   #will re-create the input file
        os.remove("SortFile.log")
        control_option = "-control " + control_file_fp
        
    shutil.copyfile(data_path,"input_data_fp.bed")
    data_path = "input_data_fp.bed"
    sort_input_file = language + space + pfms_path +"SortFiles.jar bed ./" + space + data_path
    gunzip_sorted_file = "gunzip -d " + data_path + "*.gz"
    os.system(sort_input_file)
    os.remove(data_path)
    os.system(gunzip_sorted_file)   #will re-create the input file
    os.remove("SortFile.log")
    
    if wiggle == "-wig":
        
        FindPeaks_call_command = language + space + pfms_path + "FindPeaks.jar"+ " -input" + space + data_path + space +"-name" + space + label + space + "-aligner bed -output ./" + space + control_option + space + findpeaks_options
        FindPeaks_wig_files = "./" + label + "_triangle_standard.wig.gz"
        FindPeaks_gunzip_command = "gunzip " + FindPeaks_wig_files
        
        os.system(FindPeaks_call_command)
        os.system(FindPeaks_gunzip_command)
        
        os.remove(label + "_triangle_standard.peaks")
        os.remove(label + ".log")
        shutil.copyfile(label + "_triangle_standard.wig", label + "_findpeaks.wig")
        os.remove("./" + label + "_triangle_standard.wig")
        FindPeaks_output_file = label + "_findpeaks.wig"
       
    else:
        FindPeaks_call_command = language + space + pfms_path + "FindPeaks.jar" + " -aligner bed -input" + space + data_path + space +"-name" + space + label + space + "-output ./"+ space + "-bedgraph" + space + control_option + space + findpeaks_options
        FindPeaks_bedgraph_files = "./" + label + "_triangle_standard.bedgraph.gz" 
        FindPeaks_gunzip_command = "gunzip " + FindPeaks_bedgraph_files
        
        os.system(FindPeaks_call_command)
        os.system(FindPeaks_gunzip_command)
        
        os.remove(label + "_triangle_standard.peaks")
        os.remove(label + ".log")
        output_file = open(label + "_findpeaks.bed", 'w')
        results_file = open(label + "_triangle_standard.bedgraph", 'r')
        line = results_file.readline().rstrip('\n')
        while len(line) > 0 and line[0:3] != "chr":
            output_file.write(line + '\n')
            line = results_file.readline().rstrip('\n')
        
        i = 1
        while line[0:3] == "chr":
            split_line = line.split('\t')
            output_file.write(split_line[0] + '\t' + split_line[1] + '\t' + split_line[2] + '\t' + "FP_peak_" + str(i) + '\t' + split_line[3] + '\n')
            i = i + 1
            line = results_file.readline().rstrip('\n')
        output_file.close()
        results_file.close()
        os.remove(label + "_triangle_standard.bedgraph")
        FindPeaks_output_file = label + "_findpeaks.bed"
    
    if os.path.exists(control_file_fp):
        os.remove(control_file_fp)
        os.system("rm " + label + "_triangle_standard_*")


    os.remove(data_path)
    print_message('FindPeaks is Exiting')
    return FindPeaks_output_file
