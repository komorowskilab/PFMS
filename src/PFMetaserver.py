#! /usr/bin/python2.6
'''
    @attention: Main module for the installed version of PFMS
    
    @author: Marcin Kruczyk
    @author: Husen Umer
    @organization: LCB centre for BioInformatics at Uppsala University
    @since: 15 Mar 2011
    @version: 1.1
    @note: This is the main module for installed version of PFMS (when -normal is not used during installation time).
    but if -normal was specified during installation time then another version of this module which is located in the parent directory is be used.
    
    '''

import os
import sys
import math

from src import Compare
from src import Compare_wiggle_rank
from src import FunctionsUtility
from src import Get_min_step
from src import VariableStep_Convertion
from src import ExecutePeakFinders
from src import Compare_normalized
from src import NormalizeBEDScores


def main():
    user_type ="super"
    functionsUtility = FunctionsUtility.FunctionsUtility()
    if (not functionsUtility.checkSysArguments(sys.argv[1:])):
        functionsUtility.usage(user_type)
        return 0
    
    feature_check = functionsUtility.checkFeatureAvailability(sys.argv[1:])
    if (not feature_check):
        return 0
    else:
        sys.argv = feature_check
    
    #The expected installation directory of the peakfinders and configuration file
    #If PeakFinders does not exist in the prefix dir then change the path to look at the installation dir
    if user_type=="normal":
        pfms_path = os.path.abspath("PeakFinders")
    else:
        pfms_path = sys.prefix +"/PeakFinders"
    
    #Reading the configuration file
    #pf_options is a dictionary to store peak-finder names as keys and their options as their values.
    pf_options = {'macs': "", 'hpeak': "", 'erange': "", 'findpeaks': "", 'sissr': "", 'seqsite': "", 'cisgenome': ""}
    pf_path = {'macs': "", 'hpeak': "", 'erange': "", 'findpeaks': "", 'sissr': "", 'seqsite': "", 'cisgenome': ""}
    pf_options, pf_path, java_max_memory_usage = functionsUtility.set_properties(pfms_path+"/pfms.conf")
    
    #File conversion functions to be called here
    #The function should take an input file name and act according to it's extension, if =aln then use aln2bed and return name of the new bed file
    
    #Reading and checking the input file
    name = sys.argv[(sys.argv.index("-o")) + 1]
    
    #Get the input file path
    input_data_file = os.path.abspath(sys.argv[(sys.argv.index("-i")) + 1])
    if (not os.path.exists(input_data_file)):
        print "Input file: ", input_data_file, " does not exist"
        return 0
    
    #convert input bam file to bed4
    input_type = ""
    if "-bam" in sys.argv:
        input_type = "-bam"
        input_data_file = functionsUtility.convert_to_bed(pfms_path, input_type, input_data_file, user_type)
    if "-sam" in sys.argv:
        input_type = "-sam"
        input_data_file = functionsUtility.convert_to_bed(pfms_path, input_type, input_data_file, user_type)
    
    elif "-bed" in sys.argv:
        input_type = "-bed"
    else:
        data_file_extention = input_data_file.rstrip('\n').split('.')[len(input_data_file.rstrip('\n').split('.')) - 1]
        if data_file_extention == 'bed':
            input_type = "-bed"
        elif data_file_extention == 'bam':
            input_type = "-bam"
            input_data_file = functionsUtility.convert_to_bed(pfms_path, input_type, input_data_file, user_type)
        elif data_file_extention == 'sam':
            input_type = "-sam"
            input_data_file = functionsUtility.convert_to_bed(pfms_path, input_type, input_data_file, user_type)
        else:
            print "Input file format (" + data_file_extention + ") not recognized!. Please specify input file format [-bed | -bam]"
            return 0
    
    
    
    #Append chr at the begining of each read that starts with digit
    functionsUtility.append_chr(input_data_file)
    
    #Checking the chromosome number
    #if "-all_chr" not in sys.argv:
    status = True
    if "-chr" in sys.argv:
        data_file, chr_number_sample, status = functionsUtility.filterChromosomeNumber(input_data_file, sys.argv[(sys.argv.index("-chr")) + 1])
    
    if not status:
        return 0
    
    #Reading and checking the control input file
    control_file = ""
    input_control_file = ""
    if "-control" in sys.argv:
        input_control_file = os.path.abspath(sys.argv[(sys.argv.index("-control")) + 1])
        if (not os.path.exists(input_control_file)):
            print "Control file: ", input_control_file, " does not exist."
            return 0
        #convert control bam file to bed here
        if input_type == "-bam" or input_type == "-sam":
            input_control_file = functionsUtility.convert_to_bed(pfms_path, input_type, input_control_file, user_type)
        functionsUtility.append_chr(input_control_file)
        #if "-all_chr" not in sys.argv:
        if "-chr" in sys.argv:
            control_file, chr_number, status = functionsUtility.filterChromosomeNumber(input_control_file, sys.argv[(sys.argv.index("-chr")) + 1])
        if not status:
            return 0
    #else:
    #   control_file, chr_number = functionsUtility.filterChromosoeNumber(input_control_file, str(chr_number_sample))

    wiggle = ""
    if "-wig" in sys.argv:
        wiggle = "-wig"

    #populate based on user_peakfinders (if given by user) otherwise default_peakfinders
    default_peakfinders =["macs", "cisgenome", "sissr"] #this set will be used only if user didn't specify any in sys.argv
    if '-wig' in sys.argv:
        default_peakfinders =["macs", "cisgenome", "findpeaks"]
    peakfinders = functionsUtility.peakfinders_to_be_executed(sys.argv, default_peakfinders)

    procentage = 100
    procentage_temp = 0
    #gets percentage value from user, for what ??
    if "-output_percentage" in sys.argv:
        procentage_temp = sys.argv[(sys.argv.index("-output_percentage")) + 1]
        if (procentage_temp.isdigit() and procentage_temp<=100 and procentage_temp>0):
            procentage = int(procentage_temp)
        else: print "-output_percentage value should be between 1 and 100, the default value(100) is used"

    normalization_type_bed = "-quantile75"
    normalization_type_wig = "-quantile75"
    normal_shift = 3 #only used with normal normalization

    if "-quantile" in sys.argv:
        quantile_temp = sys.argv[(sys.argv.index("-quantile")) + 1]
        if (quantile_temp.isdigit() and int(quantile_temp) <= 100 and int(quantile_temp) > 0):
            normalization_type_bed = "-quantile" + quantile_temp
            normalization_type_wig = "-quantile" + quantile_temp
        else:
            print "-output_percentage value should be between 1 and 100, the default value(100) is used"

    elif "-average" in sys.argv:
        normalization_type_bed = "-average"
        normalization_type_wig = "-average"

    elif "-rank" in sys.argv:
        normalization_type_bed = "-rank"
    #        normalization_type_wig = "-rank"

    elif "-rank_top" in sys.argv:
        normalization_type_bed = "-rank_top"

    elif "-normal" in sys.argv:
        normalization_type_bed = "-normal"
        if "-normal_shift" in sys.argv:
            normal_shift = int(sys.argv[(sys.argv.index("-normal_shift")) + 1])
    #        normalization_type_wig = "-normal"



    treshold_type = "-voting"
    if "-voting" in sys.argv:
        treshold_type = "-voting"
    if "-minFP" in sys.argv:
        treshold_type = "-minFP"
    if "-minFN" in sys.argv:
        treshold_type = "-minFN"


    Results_dir = "Results_" + name
    try:
        os.mkdir(Results_dir)
    except OSError:
        print "Creating the output directory failed."
        if os.access(Results_dir, os.F_OK) == True:
            print ("Directory already exists. Change label or remove directory: " + Results_dir)
            return 0
    os.chdir("./" + Results_dir)
    Results_dir_chrN = "."
    pfms_results_list = []
    if "-chr" not in sys.argv:
        #split the sample dataset based on chromosome number
        all_peakfinder_results = []
        input_files_list = functionsUtility.split_chromosomes(pfms_path, input_data_file, "splitted", "./", java_max_memory_usage)
        #split the control dataset based on chromosome number
        if "-control" in sys.argv:
            control_files_list = functionsUtility.split_chromosomes(pfms_path, input_control_file, "splitted_control", "./", java_max_memory_usage)
        created_directories = []    #the result of PFMS for each chromosome will be stored in it's own directory
        max_cpu_use = run_in_parallel()

        for i in range(len(input_files_list)):
            print "Processing: " + input_files_list[i]
            Results_dir_chrN = input_files_list[i][0:(len(input_files_list[i])-4)]
            os.mkdir(Results_dir_chrN)
            created_directories.append(Results_dir_chrN)
            
            if "-control" in sys.argv:
                if os.access(input_files_list[i], os.F_OK) and os.access(input_files_list[i][0:(len(input_files_list[i])-4)]+"_control.bed" , os.F_OK):
                    pfms_results_single_chr, peakfinder_results_single_chr = run_pfms("../"+input_files_list[i], "../"+input_files_list[i][0:(len(input_files_list[i])-4)]+"_control.bed", wiggle, name, pfms_path, procentage, peakfinders, user_type, pf_options, pf_path, java_max_memory_usage, Results_dir_chrN, normalization_type_bed, normalization_type_wig, normal_shift, treshold_type)
                    pfms_results_list.append(pfms_results_single_chr)
                    all_peakfinder_results.append(peakfinder_results_single_chr)
                else:
                    if not os.access(input_files_list[i], os.F_OK):
                        print "input file: "+ input_files_list[i]+ " is not exist"
                    if not os.access(input_files_list[i][0:(len(input_files_list[i])-4)]+"_control.bed" , os.F_OK):
                        print "input file: "+ input_files_list[i][0:(len(input_files_list[i])-4)]+"_control.bed"+ " is not exist"
            else:
                if max_cpu_use >1 : #run an instance of pfms for each chromosome in parallel (uses the same setting as peak-finder execution)
                    pfms_results_single_chr, peakfinder_results_single_chr = run_pfms("../"+input_files_list[i], "", wiggle, name, pfms_path, procentage, peakfinders, user_type, pf_options, pf_path, java_max_memory_usage, Results_dir_chrN, normalization_type_bed, normalization_type_wig, normal_shift, treshold_type)
                    pfms_results_list.append(pfms_results_single_chr)
                    all_peakfinder_results.append(peakfinder_results_single_chr)
                else:   #run sequentially
                    pfms_results_single_chr, peakfinder_results_single_chr = run_pfms("../"+input_files_list[i], "", wiggle, name, pfms_path, procentage, peakfinders, user_type, pf_options, pf_path, java_max_memory_usage, Results_dir_chrN, normalization_type_bed, normalization_type_wig, normal_shift, treshold_type)
                    pfms_results_list.append(pfms_results_single_chr)
                    all_peakfinder_results.append(peakfinder_results_single_chr)

        #removing the splitted sample files
        for i in range(len(input_files_list)):
            if os.access(input_files_list[i], os.F_OK):
                os.remove(input_files_list[i])
        #removing splitted chromosome files of control dataset
        if "-control" in sys.argv:
            for i in range(len(control_files_list)):
                if os.access(control_files_list[i], os.F_OK):
                    os.remove(control_files_list[i])

        #combine the results
        format = ".bed"
        if "-wig" in sys.argv:
            format = ".wig"
        pfms_final_result = name + "_all"+format
        all_results = open(pfms_final_result, "a")
        #if the WIG format was specified then the headers will be added to each result individually to specify their chromosome number.
        #but in case of BED format, chromosome number is exist for every single read so there is no need to specify them in the header of each result file.
        if '-wig' not in sys.argv:
            all_results.write("track name=\""+name+" PFMetaServer Results\""+" description=\""+name+" PFMetaServer Results\""+"\n")


        if len(pfms_results_list) > 0:
            
            for i in range(len(pfms_results_list)):
                print "appending content of: "+ pfms_results_list[i] + "to the final result"
                if os.access(pfms_results_list[i], os.F_OK):
                    file_content = open(pfms_results_list[i])
                    all_results.write(file_content.read())
                    file_content.close()
        all_results.close()
        #remove splitted files, meta_info and logo
        if "-store_results" in sys.argv:
            print "Merging peakfinder results"
            combine_results(all_peakfinder_results,name, format)

        if "-store_results" not in sys.argv and len(created_directories)>0: #delete this condition since we keep a single file for each result when user asked for that
            for i in range(len(created_directories)):
                if os.access(created_directories[i], os.F_OK) == True:
                    os.system("rm -r " +created_directories[i])

        os.remove("meta_info.txt")
        os.remove("SeparateReads.log")
        if os.access(pfms_final_result, os.F_OK):
            print "Check the output file, contains the identified peaks of all the chromosomes:" + pfms_final_result
        else:
            print "PFMS Failed, please check the output log to identify the error(s)"
    else:
        run_pfms(data_file, control_file, wiggle, name, pfms_path, procentage, peakfinders, user_type, pf_options, pf_path, java_max_memory_usage, Results_dir_chrN, normalization_type_bed, normalization_type_wig, normal_shift, treshold_type)

    if input_type == "-bam":
        os.remove(input_data_file)
        try:
            os.remove(input_control_file)
        except:
            pass


def combine_results(pf_results, name, format):
    
    pf_result_names = ["_macs", "_erange", "_cisgenome", "_sissr", "_findpeaks", "_hpeak", "_seqsite"]
    
    for i in range(len(pf_results)):
        for j in range(len(pf_results[i])):
            #print pf_results[i][j]
            for k in range(len(pf_result_names)):
                if pf_results[i][j].find(pf_result_names[k] + format) >= 0:
                    pf_result_single_chr = open(pf_results[i][j])
                    read_lines_pf_results = pf_result_single_chr.readlines()
                    copy_content = ""
                    if '-wig' not in sys.argv:
                        #remove file headers before merging - when BED format is used
                        if read_lines_pf_results[0].startswith("track"):
                            if read_lines_pf_results[1].startswith("variableStep"):
                                copy_content = "".join(read_lines_pf_results[2:])
                            else:
                                copy_content = "".join(read_lines_pf_results[1:])
                        elif read_lines_pf_results[0].startswith("browser"):
                            copy_content = "".join(read_lines_pf_results[2:])
                        elif read_lines_pf_results[0].startswith("# This BedGraph file was generated by FindPeaks"):
                            copy_content = "".join(read_lines_pf_results[3:])
                        elif read_lines_pf_results[0].startswith("# This wig file was"):
                            copy_content = "".join(read_lines_pf_results[4:])
                        else:
                            copy_content = "".join(read_lines_pf_results)
                    #But when WIG format is used, in order to know the chromosome number of each section the headers should not be removed
                    else:
                        copy_content = "".join(read_lines_pf_results)
                    del read_lines_pf_results
                    #add a header to each new file (combination result)
                    if not os.access(pf_result_names[k][1:] + "_results" + format, os.F_OK):
                        header = "track name=" + name+"_"+pf_result_names[k][1:] + "_results"+" description=" + name+"_"+pf_result_names[k][1:] + "_results"+"\n"
                    else:
                        header = ""
                    pf_last_result = open(pf_result_names[k][1:] + "_results" + format, "a")
                    pf_last_result.write(header+copy_content)
                    pf_last_result.close()
                    pf_result_single_chr.close()
                    del copy_content
                    break

def run_in_parallel():
    
    #Check whether multiprocessing module exists and more than a processor is availaable
    # Python 2.6+
    import_multiprocessing = False
    max_cpu_use = 6
    max_cpu_use_temp = 0
    min_cpu_to_work_in_parallel = 2
    min_cpu_to_work_in_parallel_temp = 0
    cpu_count = 1
    
    try:
        import multiprocessing
        import_multiprocessing = True
        cpu_count= multiprocessing.cpu_count()
        #Get and check max cpu use from user
        if '-max_cpu_use' in sys.argv:
            max_cpu_use_temp = sys.argv[(sys.argv.index("-max_cpu_use")) + 1]
            if (max_cpu_use_temp.isdigit() and max_cpu_use_temp>0): max_cpu_use = int(max_cpu_use_temp)
            else: print "-max_cpu_use value should be greater than 0, the default value (6) is used"
        else:
            max_cpu_use = cpu_count
        #Get and check min_cpu_to_work_in_parallel  from user
        if '-min_cpu' in sys.argv:
            min_cpu_to_work_in_parallel_temp = sys.argv[(sys.argv.index("-min_cpu")) + 1]
            if (min_cpu_to_work_in_parallel_temp.isdigit() and min_cpu_to_work_in_parallel_temp>0):
                min_cpu_to_work_in_parallel = int(min_cpu_to_work_in_parallel_temp)
            else:
                print "-min_cpu_to_work_in_parallel value should be greater than 0, the default value (2) is used"
    except (ImportError,NotImplementedError):
        pass
    #When -sequential not in sys.argv, multiprocessing module is available and the number of available processors was higher than min_cpu, run in parallel
    if '-sequential' in sys.argv or (cpu_count<min_cpu_to_work_in_parallel and "-parallel" not in sys.argv) or not import_multiprocessing:
        max_cpu_use = 1

    return max_cpu_use


def run_pfms(data_file, control_file, wiggle, name, pfms_path, procentage, peakfinders, user_type, pf_options, pf_path, java_max_memory_usage, Results_dir_chrN, normalization_type_bed, normalization_type_wig, normal_shift, treshold_type):
    
    os.chdir(Results_dir_chrN)
    #run the meta server and give the results
    functionsUtility = FunctionsUtility.FunctionsUtility()
    executePeakFinders = ExecutePeakFinders.ExecutePeakFinders()
    
    #when -sequential not in sys.argv and multiprocessing module is available and the number of available processors was higher than min_cpu then run in parallel
    max_cpu_use = run_in_parallel()
    pf_info = {'label':name, 'data_file':data_file,'control_file':control_file, 'wiggle':wiggle, 'pfms_path':pfms_path, 'java_max_memory_usage':java_max_memory_usage }
    peakfinder_results = [] #contains file names of all the peakfinder results (some of them might get ignored)
    if max_cpu_use>1 and len(peakfinders)>1:
        peakfinder_results = executePeakFinders.parallel(max_cpu_use, peakfinders, pf_path, pf_options, pf_info, peakfinder_results)
    else:
        peakfinder_results = executePeakFinders.sequential(peakfinders, pf_path, pf_options, pf_info, peakfinder_results)

    #Minimum file size (in Bytes) of a peakfinder result in order to be included in the comparisson
    min_file_size_to_accept_peakfinder_results = 1024 #Bytes
    if '-min_size' in sys.argv:
        min_size_temp = sys.argv[(sys.argv.index("-min_size")) + 1]
        if (min_size_temp.isdigit() and min_size_temp>=0):
            min_file_size_to_accept_peakfinder_results = int(min_size_temp)*1024    #To get it in KB
        else: print "-min_size value should be an integer, the default value (1) KB is used"

    Results = []    #contains file names of the peakfinder results that will be used in the comparison
    if len(peakfinder_results) > 0:
        for i in range(len(peakfinder_results)):
            if os.access(peakfinder_results[i], os.F_OK) and os.stat(peakfinder_results[i]).st_size >= min_file_size_to_accept_peakfinder_results:
                print peakfinder_results[i] + " : is added to the comparison"
                Results.append(peakfinder_results[i])
            else:
                print peakfinder_results[i]+" : is ignored, very few reads"

    #print "%d Peak finders will be used in the comparison" % len(Results)
    #min_rank = number of peakfinders should vote for a peak in order to be accepted as a peak
    min_rank = math.floor((len(peakfinders)/2.0)-0.1) + 1
    min_rank_temp = 0
    if '-min_rank' in sys.argv:
        min_rank_temp = sys.argv[(sys.argv.index("-min_rank")) + 1]
        
        if min_rank_temp.isdigit():
            if (int(min_rank_temp)>=0 and int(min_rank_temp) <=len(peakfinders)):
                min_rank = int(min_rank_temp)
            else:
                print "-min_rank value should be greater than 0 and smaller or equal to number of requested peak-finders, the default value (%d)is used" % min_rank
        else:
            print "-min_rank value should be a digit"
            functionsUtility.usage(user_type)
            return 0

    # if wiggle option was set then convert the results to wiggle file format
    # then compare the wiggle files otherwise compare the results directly without conversion
    output_file =  name + "_Results.bed"
    Converted_Results = []
    if len(Results) > 0:
        if wiggle == "-wig":
            step = Get_min_step.Get_min_step(Results)
            Converted_Results = []
            print "Converting the results to Wiggle VariableStep"
            for i in range(0, len(Results)):
                path, status = VariableStep_Convertion.Convert_to_VariableStep(Results[i], name, step)
                Converted_Results.append(path)
            #if Results[i] != Converted_Results[len(Converted_Results) - 1]:
            #os.remove(Results[i])
            
            print "Comparing the results..."
            Compare_wiggle_rank.Compare_results(name, Converted_Results, procentage, normalization_type_wig)
            for i in range(0, len(Converted_Results)):
                if os.access(Converted_Results[i], os.F_OK):
                    os.remove(Converted_Results[i])
        
            del Converted_Results
            output_file =  name + "_Results.wig"
        else:
            print "Comparing the results..."
            if treshold_type == "-voting":
                Compare.Compare_results(name, Results, min_rank)
            else:
                Results_normalized, none_files_with_scores = NormalizeBEDScores.NormalizeBEDScores(Results, normalization_type_bed, normal_shift)
                if none_files_with_scores:
                    print "No scores in the BED files. Starting voting procedure!"
                    Compare.Compare_results(name, Results, min_rank)
                else:
                    Compare_normalized.CompareNormalizedResults(name, Results_normalized, min_rank, treshold_type)


    if "-store_results" not in sys.argv:
        print "Removing the peak-finder results"
        for i in range(len(peakfinder_results)):
            if os.access(peakfinder_results[i], os.F_OK):
                os.remove(peakfinder_results[i])
                print "removed: " + peakfinder_results[i]
    #del Results
    del peakfinder_results

    if os.path.exists(output_file): # and os.stat(output_file).st_size !=0:
        print "PFMS executed successfully, check the output file: " + Results_dir_chrN +"/"+output_file
    else:
        print "PFMS failed, please make sure:\n- The data (and/or control) file is in the standard BED format\n- PeakFinders directory exists in the installation directory\n- Or None of the peak-finders could detect any peaks (add -store_results option to see the output of each peak-finder)"
    os.chdir("../")
    if (os.path.exists(data_file)):
        os.remove(data_file)
    if (os.path.exists(control_file)):
        os.remove(control_file)
    
    for i in range(len(Results)):
        Results[i] = Results_dir_chrN + "/" + Results[i]
    
    return Results_dir_chrN + "/" + output_file, Results

if __name__ == "__main__":
    main()
