'''

@attention: A set of utility functions, used elsewhere in this package
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 31, Jan 2011
@version: 1.0

'''

import sys
import os

class FunctionsUtility:
    '''
    This class implements the common function that will be used by other classes
    '''

    def file_exists(self,file_name):
        return os.path.exists(file_name)
    
    #check length of user input
    def checkSysArguments(self, sys_arg_list):

        PFMS_accepted_args = ["-i", "-o", "-parallel", "-sequential", "-control", "-wig", "-chr", "-help", "-bed", "-sam" , "-bam", 
                        "-macs", "-sissr", "-seqsite", "-erange", "-findpeaks", "-cisgenome", "-hpeak",
                        "-max_cpu_use", "-min_cpu", "-output_percentage", "-min_rank", "-min_size", "-store_results", "-all_chr",
                        "-voting", "-rank", "-rank_top", "-normal", "-normal_shift", "-quantile", "-average", "-minFP", "-minFN"]

        if len(sys_arg_list) >= 2 and "-i" in sys_arg_list and "-o" in sys_arg_list:
            for i in sys_arg_list:
                if i.startswith("-") and i not in PFMS_accepted_args:
                    print i + " is not in the PFMS options list"
                    return False
            return True
        else:
            return False
    
    def checkFeatureAvailability(self, sys_arg_list):
        #print "Please consider the following issues:\n"
        if ("-sissr" in sys_arg_list and "-wig" in sys_arg_list):
            print "-wig is not supported for SISSR"
            ask = raw_input("To continue the experiment without using SISSr type y otherwise type any other key to exit\t")
            if(ask == "y"):
                sys_arg_list.remove("-sissr")
            else: return False
        if ("-erange" in sys_arg_list and "-wig" in sys_arg_list):
            print "\n-wig is not supported for Erange"
            ask = raw_input("To continue the experiment without using Erange type y otherwise type any other key to exit\t")
            if(ask == "y"):
                sys_arg_list.remove("-erange")
            else: return False
        if ("-findpeaks" in sys_arg_list and not "-wig" in sys_arg_list):
            print "Findpeaks only returns bedgraph, to get fewer regions you can repeat the experiment without Findpeaks"
        if ("-output_percentage" in sys_arg_list and not "-wig" in sys_arg_list):
            print "output_percentage can only be used with -wig option"
            return False
        if ("-all_chr" in sys_arg_list and "-chr" in sys_arg_list):
            print "Please specify either -chr <number> or -all_chr"
            return False
        if ("-seqsite" in sys_arg_list and "-wig" in sys_arg_list and "-chr" not in sys_arg_list):
            print "\nSeqsite returns a very big wig file which dramatically increases the running time ."
            ask = raw_input("If you still want to use SeqSite type y otherwise type any key to continue without using SeqSite:\t")
            if(ask != "y"):
                sys_arg_list.remove("-seqsite")
        if ("-hpeak" in sys_arg_list and "-chr" in sys_arg_list):
            print "\nHPeak increases the running time dramatically."
            ask = raw_input("If you still want to use HPeak type y, otherwise press any key to continue without using HPeak:\t")
            if(ask != "y"):
                sys_arg_list.remove("-hpeak")
        if (("-bed" in sys.argv and ("-bam" in sys.argv or "-sam" in sys.argv)) or ("-bam" in sys.argv and ("-bed" in sys.argv or "-sam" in sys.argv)) or ("-sam" in sys.argv and ("-bam" in sys.argv or "-bed" in sys.argv)) ):
            print "Specify only one input type format"
            return False
        if ("-sequential" in sys_arg_list and "-parallel" in sys_arg_list):
            print "-parallel and -sequential options can't be used together"
            return False
        if ("-help" in sys_arg_list):
            return False
        if ("-wig" in sys_arg_list and ("-minFP" in sys_arg_list or "minFN" in sys_arg_list or "-voting" in sys_arg_list or "-normal" in sys_arg_list or "-rank" in sys_arg_list or "-rank_top" in sys_arg_list)):
            print "-minFP, -minFN, -voting, -normal, -rank and -rank_top are not supported for WIG"
            return False
        if (("-voting" in sys_arg_list and ("-minFP" in sys_arg_list or "minFN" in sys_arg_list)) or ("-minFN" in sys_arg_list and ("-voting" in sys_arg_list or "-minFP" in sys_arg_list)) or ("-minFP" in sys_arg_list and ("-voting" in sys_arg_list or "-minFN" in sys_arg_list))):
            print "please choose either one of these three options -minFP, -minFN or -voting"
            return False
        
        if ("-min_rank" in sys.argv):
            if sys.argv[(sys.argv.index("-min_rank")) + 1].isdigit():
                if (int(sys.argv[(sys.argv.index("-min_rank")) + 1])<0):
                    print "-min_rank value should be greater than 0 and smaller or equal to number of requested peak-finders"
                    return False
            else:
                print "-min_rank value should be a digit"
                return False
        
        if ("-normal_shift" in sys.argv):
            if ("-normal" in sys.argv):
                if sys.argv[(sys.argv.index("-normal_shift")) + 1].isdigit():
                    if (int(sys.argv[(sys.argv.index("-normal_shift")) + 1])<0):
                        print "-normal_shift value should be greater than 0 and preferably greater than 3"
                        return False
                else:
                    print "-normal_shift value should be a digit"
                    return False
            else:
                print "normal_shift can only be used with normal normalization"
                return False
        
        return sys_arg_list
        

    def usage(self,user_type):
        run_path="PFMetaserver"
        if user_type == "normal":run_path="python PFMetaserver.py"
        sys.stderr.write("\nPFMS Usage: " +
        "<"+run_path+" -i <input_file_path>  -o <output_label> [-control <control_file_path>] [-chr <chromosome>] [min_rank <number>] [-wig] [-output_percentage <number>] [-parallel | -sequential] [-max_cpu_use <number>] [-min_cpu <number>] [-voting | -minFP | minFN] [-quantile <number> | -normal | -rank | -rank_top | -average] [-macs] [-sissr] [-seqsite] [-erange] [-findpeaks] [-cisgenome] [-hpeak]>\n" +
        "Description:\n--------------------------------------------------\n" +
        "-i <input_file_path>\t\tInput data file path\n" +
        "-o <output_label>\t\t\tOutput directory and file name\n" +
        "[-control <control_file_path>]\tControl data file path\n" +
        "[-chr <chromosome>]\t\t\tTo be set when the input file includes more than one chromosome number, default(chromosome number of the first read)\n" +
        "[-all_chr]\t\t\t\tPFMS analyzes each chromosome in the given dataset individually and combines results of all the chromosomes" +
        "\n[-min_rank <number>]\t\t\tA peak is significant if it's detected by given <number> peak-finders (default is more than half)" +
        "\n[-wig] \t\t\t\t\tGives the detected peaks in WIG format while the default is BED. (please note this feature can't be used with sissr)" +
        "\n[-bed | -bam | -sam] \t\t\t\t\t\tInput file format (The tag file and control file have to have the same format"+
        "\n[-output_percentage <number>]\t\tThe percentage of the identified peaks to be obtained (default is 100)" +
        "\n[-parallel | -sequential]\t\tTo change the default specify one of them(default depends on number of available processors of the current system)"+
        "\n[-voting | -minFP | -minFN]\t\tThe peak selection methods when BED is used (minFP and minFN use the scores of the peaks obtained form the peakfinders)"+
        "\n[-quantile <number> | -normal [-normal_shift <number>]| -rank | -rank_top | -average]\t\t The normalization methods for normalizing scores of the peaks if BED is used"+
        "\n[-quantile <number> | -average]\t\t The normalization methods for normalizing scores of the peaks if WIG is used"+
        "\n[-max_cpu_use <number>]\t\t\tSets maximum number of processors (CPU) that can be used by PFMS (default is 6)" +
        "\n[-min_cpu <number>]\t\t\tPFMS is running in parallel if minimum number of processors (CPU) was available on the system (default is 2)" +
        "\n[-store_results]\t\t\t\tTo keep the original files generated by the peak-finders" +
        "\n[-min_size <number>]\t\t\tMinimum file size (in KB) of a peak-finder result to be included in the comparison (default is 1)" +
        "\n\nList of the peak-finders that can be used (the default: -macs -sissr -cisgenome)\n"+
        "[-macs] [-sissr] [-seqsite] [-erange] [-findpeaks] [-cisgenome] [-hpeak]\n\n")
        
    #populate a list based on user_peakfinders (if given by user) otherwise default_peakfinders
    def peakfinders_to_be_executed(self, sys_arg_list, default_peakfinders):

        user_peakfinders = []
        supported_peakfinders = ["-macs", "-sissr", "-seqsite", "-erange", "-findpeaks", "-cisgenome", "-hpeak"]
        for i in range(len(sys_arg_list)):
            if sys_arg_list[i] in supported_peakfinders:
                #remove the preceding (-) char and add the rest (ex: macs) to the list
                user_peakfinders.append(sys_arg_list[i][1:])

        if len(user_peakfinders) > 0:
            return user_peakfinders
        else:
            return default_peakfinders

    #Appends chr to chromosome number for all the reads that doesn't start with chr
    def append_chr(self, input_file):
        fileIn = open(input_file, "r+")
        lineList = fileIn.readlines()
        fileIn.seek(0)
        for line in lineList:
            if (line[0].isdigit() or line[0] =="X" or line[0] =="x" or line[0] =="Y" or line[0] =="y" or line[0] =="M" or line[0] =="m"):
                if line.startswith("23"): line = "chrX" + line[2:]
                elif line.startswith("24"): line = "chrY" + line[2:]
                else: line = "chr" + line
            fileIn.write(line)
        fileIn.flush()
        fileIn.close()
        return None

    def set_properties(self, config_file_path):
        '''Reading the configuration file to set peakfinders options'''

        property_file = open(config_file_path, "r")
        line_list = property_file.readlines()
        pf_options = {'macs': "", 'hpeak': "", 'erange': "", 'findpeaks': "", 'sissr': "", 'seqsite': "", 'cisgenome': ""}
        pf_path = {'macs': "", 'hpeak': "", 'erange': "", 'findpeaks': "", 'sissr': "", 'seqsite': "", 'cisgenome': ""}
        java_max_memory_usage = ""
        #A dictionary to store peak-finder names as keys and the options as their values.
        for line in line_list:
            if (line.startswith("MACS:")):
                pf_options['macs'] = line[len("MACS:"):].strip()
            elif (line.startswith("HPEAK:")):
                pf_options['hpeak'] = line[len("HPEAK:"):].strip()
            elif (line.startswith("FINDPEAKS:")):
                pf_options['findpeaks'] = line[len("FINDPEAKS:"):].strip()
            elif (line.startswith("ERANGE:")):
                pf_options['erange'] = line[len("ERANGE:"):].strip()
            elif (line.startswith("SISSR:")):
                pf_options['sissr'] = line[len("SISSR:"):].strip()
            elif (line.startswith("SEQSITE:")):
                pf_options['seqsite'] = line[len("SEQSITE:"):].strip()
            elif (line.startswith("CISGENOME: ")):
                pf_options['cisgenome'] = line[len("CISGENOME:"):].strip()

            elif (line.startswith("MACS-PATH:")):
                pf_path['macs'] = line[len("MACS-PATH:"):].strip()
            elif (line.startswith("HPEAK-PATH:")):
                pf_path['hpeak'] = line[len("HPEAK-PATH:"):].strip()
            elif (line.startswith("FINDPEAKS-PATH:")):
                pf_path['findpeaks'] = line[len("FINDPEAKS-PATH:"):].strip()
            elif (line.startswith("ERANGE-PATH:")):
                pf_path['erange'] = line[len("ERANGE-PATH:"):].strip()
            elif (line.startswith("SISSR-PATH:")):
                pf_path['sissr'] = line[len("SISSR-PATH:"):].strip()
            elif (line.startswith("SEQSITE-PATH:")):
                pf_path['seqsite'] = line[len("SEQSITE-PATH:"):].strip()
            elif (line.startswith("CISGENOME-PATH:")):
                pf_path['cisgenome'] = line[len("CISGENOME-PATH:"):].strip()
            elif (line.startswith("java_max_memory:")):
                java_max_memory_usage = line[len("java_max_memory:"):].strip()
                if not java_max_memory_usage.startswith("-X"):
                    java_max_memory_usage = ""
                else:
                    java_max_memory_usage = " " + java_max_memory_usage
            else:
                pass
        property_file.close()
        return pf_options, pf_path, java_max_memory_usage

    #Not implemented yet
    def convert_to_bed(self, pfms_path, input_type, input_file, user_type):
        #check the input_file extension
        if user_type == "normal":
            bamToBEDPath = pfms_path + "/BEDTools-Version-2.16.2/bin/bamToBed "
            samtoolsPath = pfms_path + "/samtools-0.1.18/samtools view -bS "
        else:
            bamToBEDPath = "bamToBed "
            samtoolsPath = "samtools view -bS "

        if input_type=="-bam": 
            print "Converting the input file to BED format..."
            command_convert_file = bamToBEDPath + "-i " + input_file + " > " + input_file + ".bed"
            os.system(command_convert_file)
            converted_file = input_file + ".bed"
        elif input_type=="-sam":
            print "Converting the input file to BED format..."
            converted_file_temp = input_file + "_temp.bam"
            command_sam_to_bam = samtoolsPath + input_file + " > " + converted_file_temp
            os.system(command_sam_to_bam)
            command_convert_file = bamToBEDPath + "-i " + converted_file_temp + " > " + input_file + ".bed"
            os.system(command_convert_file)
            os.remove(converted_file_temp)
            converted_file = input_file + ".bed"
        return converted_file


    def split_chromosomes(self,pfms_path, input_file, output_file_base_name,output_file_path, java_max_memory_usage):
        output_files =[]
        #zipped_file_path = os.path.dirname(input_file)+"/"
        separate_reads_tool = "java" + java_max_memory_usage + " -jar "+ pfms_path + "/SeparateReads.jar"
        print "Separating reads of the input file based on chromosome number"
        separate_command = separate_reads_tool + " bed " + input_file + " " + " " + output_file_path + " " + output_file_base_name
        os.system(separate_command)  #a zipped file for each chromosome created in Results_label folder
        #read meta data file to get the new file names and check their number of reads
        read_metadata_file = open("meta_info.txt")
        line_list = read_metadata_file.readlines()
        for line in line_list[1:]:#the first line is info
            splitted_line = line.split("\t")#first column is file name and second is number of reads
            if int(splitted_line[1]) >=10:
                zipped_file_name = splitted_line[0].rstrip()
                if os.access(zipped_file_name, os.F_OK):
                    os.system("gunzip "+output_file_path+zipped_file_name)
                    output_files.append(zipped_file_name[0:(len(zipped_file_name)-3)])
            else:
                print "Ignored: " + splitted_line[0].rstrip() + ", to few reads"
                if os.access(splitted_line[0].rstrip(), os.F_OK):
                    os.remove(splitted_line[0].rstrip())

        return output_files


    #filters the input file based on user's input if the given chr_number exists in the file otherwise based on the chr_number of the first read in the file.
    def filterChromosomeNumber(self, fileInName, chrN="-1"):
        #if chrN != x,X,y,Y,1-23 then make chrN as chrN of first read in the file
        #else filter based on value of chrN
        #stores the input file name
        status = True
        fileOutName = ""
        newFileName= os.path.basename(fileInName)#fileInName[0:fileInName.find(".bed")]
        file_path = os.path.dirname(fileInName)
        #open filein and copy all the reads for chrn into fileout
        fileIn = open(fileInName, "r")
        
        lineList = fileIn.readlines()
        
        chr = ""
        chr_number = ""
        #chrNN = "%d" % chrN
        #check if chrN is exist in the file, if yes then put value of chrN in chr and break
        for line in lineList:
            if (line.find(chrN, 3, 3+len(chrN)) != -1):
                chr = chrN
                break
            
        if chr != chrN:
            #set chr number of first read to chr
            for line in lineList:
                if line.startswith("chr", 0, 3): #read the first col
                    splited_line = line.split("\t")
                    first_column = splited_line[0].rstrip()
    #                if (first_column[3].isdigit() or first_column[3]=="X" or first_column[3]=="x" or first_column[3]=="Y" or first_column[3]=="y" or first_column[3]=="M" or first_column[3]=="m"):
                    number_of_digits = len(first_column) - len("chr")
                    index = 0
                    while index < number_of_digits:
                        chr_number += first_column[index+len("chr")]
                        index+=1
                    break
            ask = raw_input("Chromosome: "+ chrN + " does not exist in the input file, if you want to continue with: chr" + chr_number + " type y to continue otherwise press any other key to exit: ")
            if ask !='y':
                status = False 
                return fileOutName,chr_number, status         
            chr = chr_number
        
        print "Chromosome: ", chr , " is used."
        
        fileOutName =  file_path+"/"+chr +"_"+ newFileName
        fileOut = open(fileOutName, "w")
        chr = "chr" + chr
        
        for line in lineList:
            if line.startswith("chr", 0):
                splited_line = line.split("\t")
                if ((splited_line[0].rstrip()) == chr):
                    fileOut.write(line)
                    
        fileIn.close()
        fileOut.close()
        print "The input file has been checked"
        return fileOutName,chr_number, status

#Used by cisgenome
def get_cutoff_value(filename):
    #Get No_of_reads/window of the first line that has negbinomial_exp/obs < 0.1 
    c = 10
    fileIn = open(filename, "r")
    temp_line = fileIn.readline()
    while temp_line != "":
        if(temp_line[0].isdigit()):
            splited_line = temp_line.split("\t")
            if(float(splited_line[6]) < 0.1):
                c = int(splited_line[0])
                break
        temp_line = fileIn.readline()
    return c
#changes scores from float to int and adds step size to the given file
def conver_float_to_int(input_file_name, output_file_name):
    
    filein = open(input_file_name, "r")
    fileout = open(output_file_name, "w")
    linelist = filein.readlines()
    #find step size and add it to the header
    step = 25 #default step size
    i = 0
    while i <= len(linelist):
        if (linelist[i][0].isdigit() and linelist[i+1][0].isdigit()):
            linei0 = (linelist[i].split("\t"))[0].strip()   #get int value of read i and read i+1 positions
            linej0 = (linelist[i+1].split("\t"))[0].strip()
            step = int(linej0) - int(linei0)    #count step size based on the pos of read i and read i+1
            break
        i+=1
     
    #while line != "" :
    for line in linelist:
        if (line.startswith("track")):
            fileout.write(line)
        elif (line.startswith("variableStep")):
            line2 = line.rstrip("\n") + " span=%d" % step + "\n"
            fileout.write(line2)
        elif (line[0].isdigit()):
            #un-comment these lines to convert the float scores to int
            splited_line = line.split("\t")
            splited_line[1] = int(float(splited_line[1].rstrip("\n")))
            line =  "\t".join(["%s" % el for el in splited_line]) + "\n"
            fileout.write(line)

#Sort bed files by chrN then start position
def sort_bed_file(input_file, sorted_file):
    os.system("sort -k 1,1 -k2,2 -n -o "+ sorted_file + " "+ input_file)
    return sorted_file

#print information messages before and after running each peak-finder
def print_message(msg):
    print '********************'+msg+'********************'
