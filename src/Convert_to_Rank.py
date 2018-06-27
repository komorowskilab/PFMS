'''
Created on May 25, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.1

'''

import math

#convert values to ranks
def Fill_PF_File(peak_table, PF_file, normalization_value, step, max_number_of_peaks):
    PF_file_id = open(PF_file, "r")
    
    Output_file_path = PF_file[0: len(PF_file) - 4] + "_ranked" + PF_file[len(PF_file) - 4: len(PF_file)]
    Output_file_path_id = open(Output_file_path, "w")
    line = PF_file_id.readline()
     
    while line[0:1].isdigit() == False and line != "":
        Output_file_path_id.write(line)
        line = PF_file_id.readline()
        
    operating_start = 0
    last_read = 0
    operating_score = 0
    i = 0
    
    while line != "":
        data = line.split("\t")
        operating_start = int(data[0])
        operating_score = int(round(float(data[1].rstrip())))
        
        if operating_start > (last_read + step):
            while peak_table[i][0] != operating_start:
                i = i + 1
        
        rank_of_column  = int(round(float(max_number_of_peaks) * float(operating_score) / float(normalization_value)))
        Output_file_path_id.write(str(operating_start) + "\t" + str(rank_of_column) + "\n")
        last_read = operating_start
        line = PF_file_id.readline()
        
    PF_file_id.close()
    Output_file_path_id.close()
    
    return Output_file_path  


#The normal normalization and rank normalization is not correct for the wig file
'''

def Fill_PF_FileRank(peak_table, rank_table, PF_file, step):
    PF_file_id = open(PF_file, "r")
    
    Output_file_path = PF_file[0: len(PF_file) - 4] + "_ranked" + PF_file[len(PF_file) - 4: len(PF_file)]
    Output_file_path_id = open(Output_file_path, "w")
    line = PF_file_id.readline()
     
    while line[0:1].isdigit() == False and line != "":
        Output_file_path_id.write(line)
        line = PF_file_id.readline()
        
    operating_start = 0
    last_read = 0
    rank_of_peak = 0
    i = 0
    
    while line != "":
        data = line.split("\t")
        operating_start = int(data[0])
        
        if operating_start > (last_read + step):
            while peak_table[i][0] != operating_start:
                i = i + 1
                
            rank_of_peak = rank_table[peak_table[i][1]]
        
        Output_file_path_id.write(str(operating_start) + "\t" + str(rank_of_peak) + "\n")
        last_read = operating_start
        line = PF_file_id.readline()
        
    PF_file_id.close()
    Output_file_path_id.close()
    
    return Output_file_path  



def Fill_PF_FileNormal(peak_table, PF_file, average, sd, step, max_number_of_peaks):
    PF_file_id = open(PF_file, "r")
    
    Output_file_path = PF_file[0: len(PF_file) - 4] + "_ranked" + PF_file[len(PF_file) - 4: len(PF_file)]
    Output_file_path_id = open(Output_file_path, "w")
    line = PF_file_id.readline()
     
    while line[0:1].isdigit() == False and line != "":
        Output_file_path_id.write(line)
        line = PF_file_id.readline()
        
    operating_start = 0
    last_read = 0
    operating_score = 0
    i = 0
    
    
    while line != "":
        data = line.split("\t")
        operating_start = int(data[0])
        operating_score = int(round(float(data[1].rstrip())))
        
        
        if operating_start > (last_read + step):
            while peak_table[i][0] != operating_start:
                i = i + 1
        
        if sd == 0.0:
            rank_of_column = 3
        else:
            rank_of_column  = int(round(float(max_number_of_peaks) * (((float(operating_score) - float(average)) / float(sd)) + 3)))
        if rank_of_column < 0:
            rank_of_column = 0
        
        print [data[1].rstrip(), operating_score, average, sd, rank_of_column]
        Output_file_path_id.write(str(operating_start) + "\t" + str(rank_of_column) + "\n")
        last_read = operating_start
        line = PF_file_id.readline()
        
    PF_file_id.close()
    Output_file_path_id.close()
    
    return Output_file_path  



def Fill_Peaks_Rank_Table(Scores_Table, number_of_peaks, max_number_of_peaks):
    
    last_count = 0
    for i in range(0, len(Scores_Table)):
        operating_field = last_count + Scores_Table[i]
# we want to calculate the mean score for 
        Scores_Table[i] = int(round((float(max_number_of_peaks) * math.ceil((last_count + operating_field) / 2) / float(number_of_peaks)))
        last_count = operating_field
    return Scores_Table
'''

def Count_Scores(Peak_Table, max_value):

    scores = []
    for i in range(0, max_value + 1):
        scores.append(0)
    for i in range(0, len(Peak_Table)):
        scores[Peak_Table[i][1]] = scores[Peak_Table[i][1]] + 1
    return scores

def Find_max_value_in_table(table, column_number):

    max_value = 0
    for i in range(0, len(table)):
        if table[i][column_number] > max_value:
            max_value = table[i][column_number]
    return max_value


def Scan_for_peak_scores(file, step):
    
    file_id = open(file, "r")
    line = file_id.readline()
     
    while line[0:1].isdigit() == False and line != "":
        line = file_id.readline()

    Rank_Table = []    
    data = line.split("\t")
    place = int(data[0])
    score = int(round(float(data[1].rstrip())))
    sum_of_scores = 0
    number_of_peaks = 1
    current_peak_value = score
    last_place = place
    peak_start = place  

    
    while line != "":
        
        
        if place > (last_place + step):
            Rank_Table.append([peak_start, current_peak_value])
            current_peak_value = score
            peak_start = place
            number_of_peaks = number_of_peaks + 1
        
        elif score > current_peak_value:
            current_peak_value = score
            
        sum_of_scores = sum_of_scores + score
        last_place = place
        line = file_id.readline()
        if line == "":
            Rank_Table.append([peak_start, current_peak_value])
        else:
            data = line.split("\t")
            place = int(data[0])
            score = int(round(float(data[1].rstrip())))
            
    file_id.close()
    del line
    del data
    
    return Rank_Table, number_of_peaks, sum_of_scores



#This function evaluates the value of given quantile in the data
def EvaluateQuantile(scores, org_quantile_num, number_of_peaks):
    quantile_num = float(org_quantile_num) / float(100)
    previous_number_of_peaks = 0
    current_number_of_peaks = 0
    quantile_value = 0
    for i in range(0, len(scores)):
        current_number_of_peaks = current_number_of_peaks + scores[i]
        if quantile_num > float(previous_number_of_peaks) / float(number_of_peaks) and quantile_num <= float(current_number_of_peaks) / float(number_of_peaks):
            quantile_value = i
            break

    return quantile_value
    
    
def CalculateSD(scores, org_average, number_of_scores):
    sum_of_squares = float(0)
    average = float(org_average)
    for i in range(0, len(scores)):
        difference = float(scores[i]) - average
        square = difference * difference
        sum_of_squares = sum_of_squares + square
    
    sd = math.sqrt(sum_of_squares / float(number_of_scores - 1))
    return sd



def GetAllScores(scoring_table):
    all_scores = []
    for i in range(0, len(scoring_table)):
        all_scores.append(scoring_table[i][1])
        
    return all_scores


def Create_Rank_File(input_file, step, normalization_type, max_number_of_peaks):
    
#Funkcja zwraca tablice z poczatkami plikow i ich scoringami oraz liczbe peakow
    scoring_table, number_of_peaks, sum_of_scores = Scan_for_peak_scores(input_file, step)


#funkcja zwraca najwieksza wartosc w danej tablicy w danej kolumnie.
#W tym wypadku zwraca maxymalny scoring peaku

    max_value = Find_max_value_in_table(scoring_table, 1)


#Funkcja zwraca tablice w ktorej indeks to wysokosc peaku, a odpowiadajaca mu wartosc to liczba peakow o takiej 
#wysokosci w danym pliku    
    scores  = Count_Scores(scoring_table, max_value)


    '''
    if normalization_type == "-rank":
#Funkcja zamienia scoringiwyliczone przez Peak Findery na zunifikowane rankingi    
        Ranks = Fill_Peaks_Rank_Table(scores, number_of_peaks, max_number_of_peaks)
        
#Funkcja ktora wypelnia wynikowy plik z Peak Findera. header + starts + ranks    
        output_file = Fill_PF_FileRank(scoring_table, Ranks, input_file, step)
        del Ranks
#The normal normalization is not correct for the wig file
               
        elif normalization_type == "-normal":
            all_scores = GetAllScores(scoring_table)
            average = float(sum_of_scores) / number_of_peaks
            sd = CalculateSD(all_scores, average, number_of_peaks)

            output_file = Fill_PF_FileNormal(scoring_table, input_file, average, sd, step, max_number_of_peaks)
    
    else:
    '''
    if normalization_type == "-average":
         normalization_value = float(sum_of_scores / number_of_peaks)
    elif normalization_type[:9] == "-quantile":
        normalization_value =  EvaluateQuantile(scores, int(normalization_type[9:11]), number_of_peaks)

    output_file = Fill_PF_File(scoring_table, input_file, normalization_value, step, max_number_of_peaks)
    
    del scoring_table
    del scores
    
    return output_file

