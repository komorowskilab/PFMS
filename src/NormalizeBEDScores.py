'''
Created on Oct 12, 2012

@author: kruczyk
'''

import math

def GetScores(input_file):
    
    infile = open(input_file, 'r')
    
    scores = []
    no_scores = False
    line = infile.readline()
    
    header_line = 0
    while line[0:3] != "Chr" and line[0:3] != "chr" and line[0:1] != 'm' and line[0:1] != 'M' and line[0:1] != 'x' and line[0:1] != 'X' and line[0:1] != 'y' and line[0:1] != 'Y' and header_line < 20:
        if line[0:2].isdigit():
            if int(line[0:2]) > 0 and int(line[0:2]) < 23:
                break
        elif line[0:1].isdigit():
            if line[0:1] != '0':
                break
        else:
            header_line = header_line + 1
            line = infile.readline()
    
    file_contains_data = True
    if header_line < 20:             
        split_line = line.split('\t')
        
        if len(split_line) < 5:
            no_scores = True
            
            
        while line != "":
            split_line = line.split('\t')
            if no_scores:
                scores.append(1.0)
            else:
                scores.append(float(split_line[4].lstrip().rstrip()))
            line = infile.readline()
    else:
        file_contains_data = False
    
    
    infile.close()    
    
    return scores, no_scores, file_contains_data

def CalculateAverage(scores):
    sum_of_scores = 0
    number_of_scores = len(scores)
    
    for i in range(0, len(scores)):
        sum_of_scores = sum_of_scores + scores[i]
    
    average = float(sum_of_scores) / float(number_of_scores)
    
    return average

def CalculateSD(scores, average):
    sum_of_squares = 0
    number_of_scores = len(scores)
    
    for i in range(0, len(scores)):
        sum_of_squares = sum_of_squares + ((float(scores[i]) - float(average)) * (float(scores[i]) - float(average)))
    
    sd = math.sqrt(float(sum_of_squares) / float(number_of_scores - 1))
    return sd


def NormalizationNormal(input_file, average, sd, normal_shift, output_file, max_number_of_scores):
    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')
    
    
    line = infile.readline().rstrip('\n')
    
    header_line = 0
    while line[0:3] != "Chr" and line[0:3] != "chr" and line[0:1] != 'm' and line[0:1] != 'M' and line[0:1] != 'x' and line[0:1] != 'X' and line[0:1] != 'y' and line[0:1] != 'Y' and header_line < 20:
        if line[0:2].isdigit():
            if int(line[0:2]) > 0 and int(line[0:2]) < 23:
                break
        elif line[0:1].isdigit():
            if line[0:1] != '0':
                break
        else:
            header_line = header_line + 1
            line = infile.readline().rstrip('\n')
    
    file_contains_data = True
    if header_line < 20:
        while line != "":
            split_line = line.split('\t')
            if sd == 0.0:
                split_line[4] = normal_shift
            else:
                split_line[4]  = int(round(float(max_number_of_scores) * (((float(split_line[4].lstrip().rstrip()) - float(average)) / float(sd)) + normal_shift)))
            if split_line[4] < 0:
                split_line[4] = 0
            if split_line[0][0:3] == "Chr":
                split_line[0][0:3] = "chr"
            elif split_line[0][0:2].isdigit():
                split_line[0] = "chr" + split_line[0][0:2]
            elif split_line[0][0:1].isdigit():
                split_line[0] = "chr" + split_line[0][0:1]
            elif split_line[0][0:1] == 'm' or split_line[0][0:1] == 'M':
                split_line[0] = "chrM"
            elif split_line[0][0:1] == 'x' or split_line[0][0:1] == 'X':
                split_line[0] = "chrX"                
            elif split_line[0][0:1] == 'y' or split_line[0][0:1] == 'Y':
                split_line[0] = "chrY"  
            
            outfile.write(split_line[0])
            for i in range(1, len(split_line)):
                outfile.write('\t' + str(split_line[i]))
            outfile.write('\n')
            
            line = infile.readline().rstrip('\n')       
    else:
        file_contains_data = False

    infile.close()
    outfile.close()
            
    return output_file, file_contains_data
   
   
   
   
def NormalizationRank(input_file, score_clusters_ranks, output_file):
    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')
    
    line = infile.readline().rstrip('\n')
    
    header_line = 0
    while line[0:3] != "Chr" and line[0:3] != "chr" and line[0:1] != 'm' and line[0:1] != 'M' and line[0:1] != 'x' and line[0:1] != 'X' and line[0:1] != 'y' and line[0:1] != 'Y' and header_line < 20:
        if line[0:2].isdigit():
            if int(line[0:2]) > 0 and int(line[0:2]) < 23:
                break
        elif line[0:1].isdigit():
            if line[0:1] != '0':
                break
        else:
            header_line = header_line + 1
            line = infile.readline().rstrip('\n')
    
    file_contains_data = True
    if header_line < 20:
        while line != "":
            split_line = line.split('\t')
            for i in range(0, len(score_clusters_ranks)):
                if float(split_line[4].lstrip().rstrip()) == score_clusters_ranks[i][0]:
                    split_line[4] = str(score_clusters_ranks[i][1])
            if split_line[0][0:3] == "Chr":
                split_line[0][0:3] = "chr"
            elif split_line[0][0:2].isdigit():
                split_line[0] = "chr" + split_line[0][0:2]
            elif split_line[0][0:1].isdigit():
                split_line[0] = "chr" + split_line[0][0:1]
            elif split_line[0][0:1] == 'm' or split_line[0][0:1] == 'M':
                split_line[0] = "chrM"
            elif split_line[0][0:1] == 'x' or split_line[0][0:1] == 'X':
                split_line[0] = "chrX"                
            elif split_line[0][0:1] == 'y' or split_line[0][0:1] == 'Y':
                split_line[0] = "chrY"  
            
            outfile.write(split_line[0])
            for i in range(1, len(split_line)):
                outfile.write('\t' + split_line[i])
            outfile.write('\n')
            
            line = infile.readline().rstrip('\n')        
    else:
        file_contains_data = False

    infile.close()
    outfile.close()
            
    return output_file, file_contains_data

   
   
def Normalization(input_file, normalization_value, output_file, max_number_of_scores):
    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')
    
    line = infile.readline().rstrip('\n')
    
    header_line = 0
    while line[0:3] != "Chr" and line[0:3] != "chr" and line[0:1] != 'm' and line[0:1] != 'M' and line[0:1] != 'x' and line[0:1] != 'X' and line[0:1] != 'y' and line[0:1] != 'Y' and header_line < 20:
        if line[0:2].isdigit():
            if int(line[0:2]) > 0 and int(line[0:2]) < 23:
                break
        elif line[0:1].isdigit():
            if line[0:1] != '0':
                break
        else:
            header_line = header_line + 1
            line = infile.readline().rstrip('\n')
    
    file_contains_data = True
    if header_line < 20:
        while line != "":
            split_line = line.split('\t')
            split_line[4] = str(round(float(max_number_of_scores) * float(split_line[4].lstrip().rstrip()) / float(normalization_value), 2))
            if split_line[0][0:3] == "Chr":
                split_line[0][0:3] = "chr"
            elif split_line[0][0:2].isdigit():
                split_line[0] = "chr" + split_line[0][0:2]
            elif split_line[0][0:1].isdigit():
                split_line[0] = "chr" + split_line[0][0:1]
            elif split_line[0][0:1] == 'm' or split_line[0][0:1] == 'M':
                split_line[0] = "chrM"
            elif split_line[0][0:1] == 'x' or split_line[0][0:1] == 'X':
                split_line[0] = "chrX"                
            elif split_line[0][0:1] == 'y' or split_line[0][0:1] == 'Y':
                split_line[0] = "chrY"  
            
            outfile.write(split_line[0])
            for i in range(1, len(split_line)):
                outfile.write('\t' + split_line[i])
            outfile.write('\n')
            
            line = infile.readline().rstrip('\n')        
              
    else:
        file_contains_data = False

    infile.close()
    outfile.close()
            
    return output_file, file_contains_data
   
   
     

def SortTable(org_scores):
    scores = []
    for i in range(0, len(org_scores)):
        scores.append(org_scores[i])
    
    for i in range(0, len(scores)):
        for j in range(i, len(scores)):
            if scores[j] < scores[i]:
                temp = scores[i]
                scores[i] = scores[j]
                scores[j] = temp
    return scores
                   

def GenerateScoreClusters(sorted_scores):
    clusters = []
    scores_values = []
    scores_peaks_numbers = []
    for i in range(0, len(sorted_scores)):
        if sorted_scores[i] in scores_values:
            scores_peaks_numbers[scores_values.index(sorted_scores[i])] = scores_peaks_numbers[scores_values.index(sorted_scores[i])] + 1
        else:
            scores_values.append(sorted_scores[i])
            scores_peaks_numbers.append(1)
    
    for i in range(0, len(scores_values)):
        clusters.append([scores_values[i], scores_peaks_numbers[i]])
    
    return clusters
    


def GenerateRanksForClusters(score_clusters, number_of_scores, max_number_of_scores):    
    last_count = 0
    for i in range(0, len(score_clusters)):
        operating_field = last_count + score_clusters[i][1]
# we want to calculate the mean score for all peaks with the same score
        score_clusters[i][1] = int(round(max_number_of_scores * float(last_count + operating_field) / 2 / float(number_of_scores)))
#        score_clusters[i][1] = int(round(((last_count + operating_field) / 2) + max_number_of_scores - number_of_scores))
        last_count = operating_field
    return score_clusters
    

def GenerateRanksForClustersTop(score_clusters, max_number_of_scores):    
    last_count = 0
    for i in range(0, len(score_clusters)):
        operating_field = last_count + score_clusters[len(score_clusters) - 1 - i][1]
# we want to calculate the mean score for all peaks with the same score
        score_clusters[len(score_clusters) - 1 - i][1] = int(round(max_number_of_scores - float(last_count + operating_field - 1) / 2))
#        score_clusters[i][1] = int(round(((last_count + operating_field) / 2) + max_number_of_scores - number_of_scores))
        last_count = operating_field
    return score_clusters


def CalculateQuantile(sorted_scores, org_quantile_num):
    quantile_num = float(org_quantile_num) / float(100)
    number_of_scores = len(sorted_scores)
    index_of_quantile = int(round(quantile_num * number_of_scores)) - 1
    quantile_value = sorted_scores[index_of_quantile]
    return quantile_value



def SetMissingScores(input_file, output_file, average_score):
    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')
    
    line = infile.readline().rstrip('\n')
    
    header_line = 0
    while line[0:3] != "Chr" and line[0:3] != "chr" and line[0:1] != 'm' and line[0:1] != 'M' and line[0:1] != 'x' and line[0:1] != 'X' and line[0:1] != 'y' and line[0:1] != 'Y' and header_line < 20:
        if line[0:2].isdigit():
            if int(line[0:2]) > 0 and int(line[0:2]) < 23:
                break
        elif line[0:1].isdigit():
            if line[0:1] != '0':
                break
        else:
            header_line = header_line + 1
            line = infile.readline().rstrip('\n')
    
    if header_line < 20:
        peak_number = 0
        while line != "":
            split_line = line.split('\t')
            peak_number = peak_number + 1
            if len(split_line) == 3:
                split_line.append("peak_" + str(peak_number))
            split_line.append(str(round(float(average_score), 2)))
            if split_line[0][0:3] == "Chr":
                split_line[0][0:3] = "chr"
            elif split_line[0][0:2].isdigit():
                split_line[0] = "chr" + split_line[0][0:2]
            elif split_line[0][0:1].isdigit():
                split_line[0] = "chr" + split_line[0][0:1]
            elif split_line[0][0:1] == 'm' or split_line[0][0:1] == 'M':
                split_line[0] = "chrM"
            elif split_line[0][0:1] == 'x' or split_line[0][0:1] == 'X':
                split_line[0] = "chrX"                
            elif split_line[0][0:1] == 'y' or split_line[0][0:1] == 'Y':
                split_line[0] = "chrY"  
            
            outfile.write(split_line[0])
            for i in range(1, len(split_line)):
                outfile.write('\t' + split_line[i])
            outfile.write('\n')
            
            line = infile.readline().rstrip('\n')      
        
    infile.close()
    outfile.close()
            
    return output_file 



def GetAverageScore(input_file):
    infile = open(input_file, 'r')
    
    line = infile.readline().rstrip('\n')
    
    sum_of_scores = 0.0
    number_of_peaks = 0
    while line != "":
        split_line = line.split('\t')
        sum_of_scores = sum_of_scores + float(split_line[4].lstrip().rstrip())
        number_of_peaks = number_of_peaks + 1
        line = infile.readline().rstrip('\n')
    
    average = sum_of_scores / number_of_peaks
    return average 



def NormalizeBEDFile(file, normalization_type, normal_shift, index, max_number_of_scores):
    scores, no_scores, file_contains_data = GetScores(file)
    index = index + 1
    output_file = ""
    if normalization_type == "-normal":
        average = CalculateAverage(scores)
        sd = CalculateSD(scores, average)
        output_file, file_contains_data = NormalizationNormal(file, average, sd, normal_shift, "BED_norm_" + str(index) + ".bed", max_number_of_scores)
    elif normalization_type == "-rank" or normalization_type == "-rank_top":
        sorted_scores = SortTable(scores)
        number_of_scores = len(scores)
        score_clusters = GenerateScoreClusters(sorted_scores)
        if normalization_type == "-rank":
            score_clusters_ranks = GenerateRanksForClusters(score_clusters, number_of_scores, max_number_of_scores)
        elif normalization_type == "-rank_top":
            score_clusters_ranks = GenerateRanksForClustersTop(score_clusters, max_number_of_scores)
        output_file, file_contains_data = NormalizationRank(file, score_clusters_ranks, "BED_norm_" + str(index) + ".bed")
    elif normalization_type == "-average":
        average = CalculateAverage(scores)
        output_file, file_contains_data = Normalization(file, average, "BED_norm_" + str(index) + ".bed", max_number_of_scores)
    elif normalization_type[0:9] == "-quantile":
        sorted_scores = SortTable(scores)
        quantile = CalculateQuantile(sorted_scores, int(normalization_type[9:11]))
        output_file, file_contains_data = Normalization(file, quantile, "BED_norm_" + str(index) + ".bed", max_number_of_scores)
    
    return output_file
                          
        
      

def NormalizeBEDScores(files, normalization_type, normal_shift):
    output_files = []
    max_number_of_scores = 0
    for i in range(0, len(files)):
        scores, no_scores, file_contains_data = GetScores(files[i])
        if max_number_of_scores < len(scores):
            max_number_of_scores = len(scores)
    
    files_without_scores = []
    none_files_with_scores = True
    for i in range(0, len(files)):
        scores, no_scores, file_contains_data = GetScores(files[i])
        if file_contains_data:
            if no_scores:
                files_without_scores.append(files[i])
            else:
                none_files_with_scores = False
                index = len(output_files)
                output_files.append(NormalizeBEDFile(files[i], normalization_type, normal_shift, index, max_number_of_scores))
    

    if none_files_with_scores:
        output_files = files
    else:
        average_score = GetAverageScore(output_files[0])

        
        for i in range(0, len(files_without_scores)):
            index = len(output_files)
            output_files.append(SetMissingScores(files_without_scores[i], "BED_norm_" + str(index + 1) + ".bed", average_score))
        
    return output_files, none_files_with_scores
                
        
        
        