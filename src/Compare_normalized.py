'''
Created on Oct 12, 2012

@author: kruczyk
'''


'''
Created on Apr 20, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.0

'''
import sys, os, NormalizeBEDScores

#Funkcja sortujaca numery oznaczajace poczatki i konce plikow.
#Zwraca tablice posortowanych numerow, 
#tablice z odpowiadajacymi im znacznikami czy to poczatek czy koniec
#oraz tablice z numerem pliku wejsciowego, z ktorego numer zostal pobrany

# sort function numbers indicate the beginnings and ends of files.
# Return an array of sorted numbers
# array with the corresponding markers either the beginning or end
# and the number plates of the input file, from which the number has been downloaded

def Sort_numbers(numbers, file_numbers, markers, scores):

    index = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    for i in range(0, len(numbers)):

        min_number = numbers[i]
        min_marker = markers[i]
        min_number_file = file_numbers[i]
        min_score = scores[i]
        index = i
        for j in range(i, len(numbers)):
            if len(min_number) > len(numbers[j]):
                min_number = numbers[j]
                min_marker = markers[j]
                min_number_file = file_numbers[j]
                min_score = scores[j]
                index = j

            elif len(min_number) == len(numbers[j]):
                if min_number > numbers[j]:
                    min_number = numbers[j]
                    min_marker = markers[j]
                    min_number_file = file_numbers[j]
                    min_score = scores[j]
                    index = j
            else:
                pass
        
        tmp1 = numbers[i]
        tmp2 = markers[i]
        tmp3 = file_numbers[i]
        tmp4 = scores[i]
        
        numbers[i] = min_number
        markers[i] = min_marker
        file_numbers[i] = min_number_file
        scores[i] = min_score
        
        numbers[index] = tmp1
        markers[index] = tmp2
        file_numbers[index] = tmp3
        scores[index] = tmp4
        
    return numbers, file_numbers, markers, scores
        

def CompareNormalizedResults(label, files, min_rank, treshold_type):
    results_files_id = []
 
 
#Otwieram pliki z danymi   

    number_of_Peak_Finder_files = len(files)
    
    for i in range(0,number_of_Peak_Finder_files):
        results_files_id.append(open(files[i], "r"))

#    EOF_markers = []
#    i = 0
#    while len(results_files_id) >= i:
#        EOF_markers.append(0)
#        i = i + 1
    
#Tworze i otwieram plik wynikowy
    output_file_name = "./" + label + "_Results.bed"
    output_file_id = open(output_file_name, "w")
    
#Tworze i otwieram plik pomocniczy (tymczasowy)
    temp_compare_file = "./temp_compare.dat"
    temp_compare_id = open(temp_compare_file, "w")

#Indeks tabelki line oznacza tylko numer pliku z ktorego zostala przeczytana linia, a nie numer linii
#Index table line is only a file number from which it was read a line, not a line number    
    line = []

#Zmienna oznaczajaca liczbe blikow wejsciowych, ktorych nie udalo sie odczytac.
#Powinna byc rowna 0.
#The variable representing the number lished the input, which could not be read.
# Should be equal to 0
    unread_files_numbers = []

#W tej petli ustawiam sie na pierwszej linii z danymi. Omijam caly naglowek
# In this loop moves to the first line of data. I get around the whole header
#removes all the lines that do not has chr at the begining
    for i in range(0,number_of_Peak_Finder_files):
        line.append(results_files_id[i].readline())
        substring = line[i][0:3]
        unrecognized_file_break = 0
        while substring != "chr" and substring != "Chr" and unrecognized_file_break < 20:
            line[i] = results_files_id[i].readline()
            substring = line[i][0:3]
            unrecognized_file_break = unrecognized_file_break + 1
            if unrecognized_file_break == 19:
                unread_files_numbers.append(i)
    
    
#W tej petli ustawiam sie na pierwszej linii z danymi, 
#ale dla plikow ktore maja inny format danych. 
#Tzn. na poczatku maja tylko numer chromosomu.

#*************************************************************************
#Trzeba byloby jeszcze sprawdzic czy zakres numerow jest od 1-23 i X i Y
#*************************************************************************
#file that couldn't be read
    finally_unread_files_numbers = []
    for i in range(0,len(unread_files_numbers)):
        results_files_id[unread_files_numbers[i]].close()
        results_files_id[unread_files_numbers[i]] = open(files[unread_files_numbers[i]], "r")
        line[unread_files_numbers[i]] = results_files_id[unread_files_numbers[i]].readline()
        substring_untypical = line[unread_files_numbers[i]][0:1]
        substring_x_or_y = line[unread_files_numbers[i]][0]
        unrecognized_read_count = 0
        while substring_untypical.isdigit() == False and substring_x_or_y.isdigit() == False and substring_x_or_y != "x" and substring_x_or_y != "X" and substring_x_or_y != "y" and substring_x_or_y != "Y" and substring_x_or_y != "m" and substring_x_or_y != "M" and unrecognized_read_count < 20:
            line[unread_files_numbers[i]] = results_files_id[unread_files_numbers[i]].readline()
            substring_untypical = line[unread_files_numbers[i]][0:2]
            substring_x_or_y = line[unread_files_numbers[i]][0:1]
            unrecognized_read_count = unrecognized_read_count + 1
            if unrecognized_read_count == 19:
                finally_unread_files_numbers.append(unread_files_numbers[i])
    del unread_files_numbers
    
#**************************************************************
#Trzeba sprawdzic czy ten sam chromosom we wszystkich plikach 
#**************************************************************

    
#Tablica w ktorej przechowywane sa numery poczatku lub konca peaku
#stores start and end of each read
    Begin_End_array = []

#Tablica w ktorej przechowywane sa informacje, czy dany numer jest poczatkiem czy koncem peaku
#start with the first 
    Begin_End_array_marker = []
    
#Tablica w ktorej przechowywana jest informacja o numerze pliku z ktorego pochodzi numer
#(poczatka lub konca peaku)
#stores for wich file is it.    
    Begin_End_file_numbers_array = []
    
    
    start = "start"
    stop = "stop"

#In this table scores of the peaks are stored
    Begin_End_scores_array = []

#W tej petli dziele linie odczytane z plikow (pierwsze linie z danymi)
#i dane wpisuje do odpowiednich tabeli (zadeklarowanych powyzej)    
#fill the first table
    for i in range(0, number_of_Peak_Finder_files):
        splited_results = line[i].split('\t')
        Begin_End_array.append(splited_results[1])
        Begin_End_array_marker.append(start)
        Begin_End_file_numbers_array.append(i)
        Begin_End_scores_array.append(float(splited_results[4]))
        Begin_End_array.append(splited_results[2])
        Begin_End_array_marker.append(stop)
        Begin_End_file_numbers_array.append(i)
        Begin_End_scores_array.append(float(splited_results[4]))


#Sprawdzam numer chromosomu i go odczytuje. Nie mozna jako pierwszego ustawic innego
#Peak Findera poniewaz moze miec inny format numeru chromosomu.    
#checking number of chromosome.
    if line[0][1] == '\t' or line[0][1] == ' ':
        chromosome = "chr" + line[0][0:1]
    elif line[0][2] == '\t' or line[0][2] == ' ':
        chromosome = "chr" + line[0][0:2]    
    elif line[0][4] == '\t' or line[0][4] == ' ':
        chromosome = "chr" + line[0][3:4]
    else:
        chromosome = "chr" + line[0][3:5]  
        
    
 
#    print(line[0][4])
    
    
#Tablica z ustawionymi markerami konca pliku na 1. Jak to zostanie przestawione na 0, 
#to plik nie bedzie juz dalej odczytywany.
#stores info about the files under analysis. if it's still not done then it's value is 1 otherwise 0
    EOF_markers = []
    
    
    for i in range(0, number_of_Peak_Finder_files):
        EOF_markers.append(1)

#Zmienna pokazujaca liczbe odczytywanych nadal plikow. Jak spadnie do 0, program sie zatrzyma.        
    End_of_Files_marker = number_of_Peak_Finder_files
    
    
    new_peak_score = float(0)
    peak_votes_number = 0
    
#Sortujemy pobrane za pierwszym razem wartosci    
#sort based on the first array and change value of other arrays accordingly
    Begin_End_array, Begin_End_file_numbers_array, Begin_End_array_marker, Begin_End_scores_array = Sort_numbers(Begin_End_array, Begin_End_file_numbers_array, Begin_End_array_marker, Begin_End_scores_array)

#comparison loop
#Ta petla wypelnia plik tymczasowy    
    while End_of_Files_marker > 0:
        
#Sortuje tablice z numerami, ze znaczikami start i stop i z numerami plikow z ktorych pochodza dane        
        
#Pobieram dane z kompletem informacji (start/stop i numer pliku) i usuwam ja od razu z tablicy        
        operation_number = Begin_End_array.pop(0)   #take the first number from first array
        operation_marker = Begin_End_array_marker.pop(0)    
        
#Zmienna przechowujaca numer pliku, ktorym sie aktualnie zajmujemy        
        operation_file_number = Begin_End_file_numbers_array.pop(0)
        operation_score = float(Begin_End_scores_array.pop(0))
        
#Jesli start, to po prostu wpisze do pliku podnoszac jednoczesnie liczbe zliczen        
#k means number of overlapped peaks
#if it's start pos the increase otherwise decrease
        if operation_marker == "start":
            new_peak_score = new_peak_score + operation_score
            peak_votes_number = peak_votes_number + 1
            
#Jesli stop, to poza zminiejszeniem liczby zliczen (o jedna pozycje za wczesnie) musze jeszcze zaczytac
#kolejna dana z pliku, z ktorego wiersz wlasnie skonczylismy
        elif operation_marker == "stop":
            peak_votes_number = peak_votes_number - 1
            if peak_votes_number == 0:
                new_peak_score = 0
            else:
                new_peak_score = new_peak_score - operation_score

            
            line[operation_file_number] = ""
            line[operation_file_number] = results_files_id[operation_file_number].readline() 
            
#Sprawdzam czy to nie byl koniec pliku
            if line[operation_file_number] != "":
                splited_results = line[operation_file_number].split("\t")
                Begin_End_array.append(splited_results[1])
                Begin_End_array_marker.append(start)
                Begin_End_file_numbers_array.append(operation_file_number)
                Begin_End_scores_array.append(float(splited_results[4]))
                Begin_End_array.append(splited_results[2])
                Begin_End_array_marker.append(stop)
                Begin_End_file_numbers_array.append(operation_file_number)
                Begin_End_scores_array.append(float(splited_results[4]))
                Begin_End_array, Begin_End_file_numbers_array, Begin_End_array_marker, Begin_End_scores_array = Sort_numbers(Begin_End_array, Begin_End_file_numbers_array,Begin_End_array_marker, Begin_End_scores_array)                
            else:
                EOF_markers[operation_file_number] = 0
                End_of_Files_marker = 0
                for i in range(0, number_of_Peak_Finder_files):
                    End_of_Files_marker = End_of_Files_marker + EOF_markers[i]

#Wpisuje dane do pliku tymczasowego                     
        temp_compare_id.write(operation_number + "\t" + operation_marker + "\t" + str(new_peak_score) + '\t' + str(peak_votes_number) + "\n")
#        print(operation_number + "\t" + operation_marker + "\t" + str(new_peak_score))


#Zamykam wszystkie pliki wejsciowe
    del line
    for i in range(0,number_of_Peak_Finder_files):              
        results_files_id[i].close()
    
    del results_files_id
    del EOF_markers
    del Begin_End_array
    del Begin_End_array_marker
    del Begin_End_file_numbers_array
    del Begin_End_scores_array
    del splited_results

#Zamykam plik tymczasowy
#tmp_file stores all the information about all the reads in different files (opeation number, operation marker, and number of overlapped reads)
    temp_compare_id.close()
    
#Otwieram plik tymczasowy do czytania   
    
    temp_compare_id = open(temp_compare_file, "r")
    if "-chr" in sys.argv:
        output_file_id.write("track name=\""+label+" PFMetaServer Results\" description=\""+label+" PFMetaServer Results\""+"\n")
    
#Odczytuje pierwsza linie z pliku tymczasowego    
    temp_line = temp_compare_id.readline().rstrip('\n')

#W tej petli ustalane sa potencjalne tresholdy najmniejszy z wystaeczajaca liczba glosow i najwiekszy bez wystarczajacej liczby glosow
    largest_without_min_rank = 0
    smallest_with_min_rank = -1.0

    while temp_line != "":
        splited_line = temp_line.strip().split("\t")
        if int(splited_line[3]) >= min_rank and (float(splited_line[2]) < smallest_with_min_rank or smallest_with_min_rank == -1.0):
            smallest_with_min_rank = float(splited_line[2])
        if int(splited_line[3]) < min_rank and float(splited_line[2]) > largest_without_min_rank:
            largest_without_min_rank = float(splited_line[2])
        temp_line = temp_compare_id.readline().rstrip('\n')


    
    if treshold_type == "-minFP":
        treshold = largest_without_min_rank
    elif treshold_type == "-minFN":
        treshold = smallest_with_min_rank



#Ustawiam sie znow na poczatku pliku
    temp_compare_id.seek(0)
#Odczytuje pierwsza linie z pliku tymczasowego    
    temp_line = temp_compare_id.readline()        

# W tej petli przeprowadzone jest glosowanie i uzupelniany jest plik wynikowy BED
#at least half of the peakfinders should has found it

    max_score_of_the_peak = 0
    peak_number = 0
    we_are_in_peak = False
    while temp_line != "":
        splited_line = temp_line.strip().split("\t")
        if splited_line[1] == "start" and splited_line[3].isdigit():
            if treshold_type == "-minFP":
                if we_are_in_peak:
                    if float(splited_line[2]) > max_score_of_the_peak:
                        max_score_of_the_peak = float(splited_line[2])                    
                if float(splited_line[2]) > treshold and not we_are_in_peak: #default value is: (math.floor(len(peakfinders)/2-0.1) + 1)
                    if float(splited_line[2]) > max_score_of_the_peak:                        
                        max_score_of_the_peak = float(splited_line[2])  
                    output_file_id.write(chromosome + "\t" + splited_line[0] + "\t")
                    we_are_in_peak = True
            elif treshold_type == "-minFN":
                if we_are_in_peak:
                    if float(splited_line[2]) > max_score_of_the_peak:
                        max_score_of_the_peak = float(splited_line[2]) 
                if float(splited_line[2]) >= treshold and treshold > -1.0 and not we_are_in_peak: #default value is: (math.floor(len(peakfinders)/2-0.1) + 1)
                    if float(splited_line[2]) > max_score_of_the_peak:
                        max_score_of_the_peak = float(splited_line[2]) 
                    output_file_id.write(chromosome + "\t" + splited_line[0] + "\t")
                    we_are_in_peak = True
                    

#            print(chromosome + "\t" + splited_line[0] + "\t")
        elif splited_line[1] == "stop":
            if float(splited_line[2]) < treshold and we_are_in_peak: #default value is: math.floor(number_of_Peak_Finder_files/2-0.1)
                peak_number = peak_number + 1    
                output_file_id.write(splited_line[0] + '\t' + "PFMS_peak_" + str(peak_number) + '\t' + str(max_score_of_the_peak) + "\n")
                max_score_of_the_peak = 0
                we_are_in_peak = False
            elif float(splited_line[2]) == treshold and we_are_in_peak and treshold_type == "-minFP": #default value is: math.floor(number_of_Peak_Finder_files/2-0.1)
                peak_number = peak_number + 1    
                output_file_id.write(splited_line[0] + '\t' + "PFMS_peak_" + str(peak_number) + '\t' + str(max_score_of_the_peak) + "\n")
                max_score_of_the_peak = 0
                we_are_in_peak = False

#            print(splited_line[0] + "\n")
        temp_line = ""
        
        temp_line = temp_compare_id.readline()
#        print temp_line     

#Zamykam plik tymczasowy i plik wynikowy    
    temp_compare_id.close()
    output_file_id.close() 
    
    del splited_line
    
    #os.remove(temp_compare_file)
    #for i in range(0, len(files)):
        #os.remove(files[i])


    return 0

    