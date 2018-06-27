'''
Created on May 21, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.0

'''
import os
from src import Header_Information
from src import Output_File_Fill
from src import Temp_File_Fill
from src import Convert_to_Rank
from src import Select_Best_Peaks

def Compare_results(label, files, procentage, normalization_type):
    
    status = 1
    step = -1
    chromosome = ""
    i = 0
    
    for i in range(0, len(files)):
        status, temp_step, temp_chromosome = Header_Information.Get_Header_Info(files[i])
        if status == 0:
            if chromosome == "" and temp_chromosome != "":
                chromosome = temp_chromosome
            if temp_step > 0:
                if step <= 0:
                    step = temp_step
                else:
                    if temp_step < step:
                        step = temp_step
    
    Ranked_files = []
    max_number_of_peaks = 0
    for i in range(0, len(files)):
        scoring_table, number_of_peaks, sum_of_scores = Convert_to_Rank.Scan_for_peak_scores(files[i], step)
        if number_of_peaks > max_number_of_peaks:
            max_number_of_peaks = number_of_peaks
        
    for i in range(0, len(files)):
        Ranked_files.append(Convert_to_Rank.Create_Rank_File(files[i], step, normalization_type, max_number_of_peaks))
    
    temp_file_path = Temp_File_Fill.Fill_Temp_File(Ranked_files, label, step)
    
    for i in range(0, len(Ranked_files)):
        if os.access(Ranked_files[i], os.F_OK):
            os.remove(Ranked_files[i])
    
    full_results_path = Output_File_Fill.Fill_Output_File(temp_file_path, label, chromosome, step)
    if os.access(temp_file_path, os.F_OK):
        os.remove(temp_file_path)

    final_peak_ranking_table, number_of_peaks, sum_of_scores = Convert_to_Rank.Scan_for_peak_scores(full_results_path, step)

    max_value = Convert_to_Rank.Find_max_value_in_table(final_peak_ranking_table, 1)

    Ranks = Convert_to_Rank.Count_Scores(final_peak_ranking_table, max_value)
    
    if procentage < 100:
        Select_Best_Peaks.Select_Highest_Peaks(final_peak_ranking_table, Ranks, full_results_path, number_of_peaks, step, procentage, label, chromosome)
    
    del Ranked_files
    del final_peak_ranking_table
    del Ranks
    
    return 0


    
        