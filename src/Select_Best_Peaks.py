'''
Created on May 28, 2010
@author: kruczyk
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 14 Jan 2011, some fixes
@version: 1.0

'''

def Find_Cut_Value(Ranks_Table, number_of_peaks, procentage):
    '''Get a subset of the best ranked peaks based on -output_percantage value given by user'''

    number_of_selected_peaks = int(round((float(number_of_peaks) * float(procentage)) / 100))
    counted_peaks = 0
    i = len(Ranks_Table) - 1
    while counted_peaks < number_of_selected_peaks and i >= 0:
        counted_peaks = counted_peaks + Ranks_Table[i]
        i = i - 1
    return i

def Fill_PF_File(Beginings_Rankings_Table, Cut_Value, results_file, step, label, procentage, chromosome):

    results_file_id = open(results_file, "r")
    output_file_name = "./" + label+ "_" + str(procentage) + "_proc_results.wig"
    output_file_id = open(output_file_name, "w")
    
    header = "track type=wiggle_0 name=\"" + label + "_" + str(procentage) + "_PFMeataserver\" description=\"Peak Finder Metaserver " + str(procentage) +"% results\" color=50,50,150 visibility=2\nvariableStep chrom=" + chromosome + " span=" + str(step) + "\n"
    output_file_id.write(header)
    line = results_file_id.readline()
    while line[0:1].isdigit() == False and line != "":
        line = results_file_id.readline()
        
    operating_start = 0
    last_read = 0
    operating_score = 0
    rank_of_peak = 0
    i = 0

    while line != "":
        data = line.split("\t")
        operating_start = int(data[0])
        operating_score = int(data[1])
        
        if operating_start > (last_read + step):
            while Beginings_Rankings_Table[i][0] != operating_start:
                i = i + 1
       
            rank_of_peak = Beginings_Rankings_Table[i][1]
        if rank_of_peak >= Cut_Value:
            output_file_id.write(str(operating_start) + "\t" + str(operating_score) + "\n")
        last_read = operating_start
        line = results_file_id.readline()
        
    results_file_id.close()
    output_file_id.close()
    
    return output_file_name  


def Select_Highest_Peaks(Beginings_Rankings_Table, Ranks_Table, results_file, number_of_peaks, step, procentage, label, chromosome):
    
    Cut_Value = Find_Cut_Value(Ranks_Table, number_of_peaks, procentage)
    Output_procentage_file = Fill_PF_File(Beginings_Rankings_Table, Cut_Value, results_file, step, label, procentage, chromosome)
    
    return Output_procentage_file
    
    
