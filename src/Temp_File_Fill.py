'''
Created on May 21, 2010

@author: kruczyk
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 25 Mar 2011, some fixes
@version: 1.0

'''


def Sort_table_by_column(table, column_number):       
            
    index = 0

    for i in range(0, len(table)):

        min_value = table[i][column_number]
        index = i
        for j in range(i, len(table)):
            if min_value > table[j][column_number]:
                min_value = table[j]
                index = j
            
        temp = table[i]
        table[i] = table[index]
        table[index] = temp
    return table
                 
#Funkcja sortujaca numery oznaczajace poczatki i konce plikow.
#Zwraca tablice posortowanych numerow, 
#tablice z odpowiadajacymi im znacznikami czy to poczatek czy koniec
#oraz tablice z numerem pliku wejsciowego, z ktorego numer zostal pobrany

def Sort_numbers(numbers, file_numbers, markers, start_stops):
    

    index = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    for i in range(0, len(numbers)):

        min_number = numbers[i]
        min_marker = markers[i]
        min_number_file = file_numbers[i]
        min_start_stop = start_stops[i]
        index = i
        for j in range(i, len(numbers)):
            if len(min_number) > len(numbers[j]):
                min_number = numbers[j]
                min_marker = markers[j]
                min_number_file = file_numbers[j]
                min_start_stop = start_stops[j]
                index = j

            elif len(min_number) == len(numbers[j]):
                if min_number > numbers[j]:
                    min_number = numbers[j]
                    min_marker = markers[j]
                    min_number_file = file_numbers[j]
                    min_start_stop = start_stops[j]
                    index = j
            else:
                pass
            
        tmp1 = numbers[i]
        tmp2 = markers[i] 
        tmp3 = file_numbers[i]
        tmp4 = start_stops[i]
        numbers[i] = min_number
        markers[i] = min_marker
        file_numbers[i] = min_number_file
        start_stops[i] = min_start_stop
        numbers[index] = tmp1
        markers[index] = tmp2
        file_numbers[index] = tmp3
        start_stops[index] = tmp4
        
    return numbers, file_numbers, markers, start_stops


def Fill_Temp_File(files, label, step):
    
    results_files_id = []
 
#Otwieram pliki z danymi   OPen data files

    number_of_Peak_Finder_files = len(files)
    
    for i in range(0,number_of_Peak_Finder_files):
        results_files_id.append(open(files[i], "r"))

#    EOF_markers = []
#    i = 0
#    while len(results_files_id) >= i:
#        EOF_markers.append(0)
#        i = i + 1
    
    
#Tworze i otwieram plik pomocniczy (tymczasowy)		create and open temp file
    temp_wig_compare_file = "./" + label + "_temp_wig_compare.dat"
    temp_wig_compare_id = open(temp_wig_compare_file, "w")

#Indeks tabelki line oznacza tylko numer pliku z ktorego zostala przeczytana linia, a nie numer linii    
# index of line table means number of file (from which file you have gotten the line, line1 is from file 1,..)
    line = []

#Zmienna oznaczajaca liczbe plikow wejsciowych, ktorych nie udalo sie odczytac.
#Powinna byc rowna 0.
#IO file error. if you couldn't read the file.
    unread_files_numbers = []


#W tej petli ustawiam sie na pierwszej linii z danymi. Omijam caly naglowek  
#omit the header of data files, then start reading

    for i in range(0, number_of_Peak_Finder_files):
        line.append(results_files_id[i].readline())
        if line == "":
            unread_files_numbers.append(i)
            line.pop(i)
            continue
        
            
        substring = line[i][0:1]
        unrecognized_file_break = 0

        while substring.isdigit() == False or line[i] == "":
         
            line[i] = results_files_id[i].readline()
            substring = line[i][0:1]
            unrecognized_file_break = unrecognized_file_break + 1
            if unrecognized_file_break == 19:
                unread_files_numbers.append(i)
                line.pop(i)

  
  
    #to do
#*************************************************************************
#Trzeba byloby jeszcze sprawdzic czy zakres numerow jest od 1-23 i X i Y
#*************************************************************************
# check range of numbers of chromosom if it is between 1-23 or x or y 
    
#**************************************************************
#Trzeba sprawdzic czy ten sam chromosom we wszystkich plikach 
#**************************************************************
#check if all the files have the same chromosom
    
#Tablica w ktorej przechowywane sa numery poczatka peaku
#store number of begining and end of peaks
    Begin_End_array = []

#Tablica w ktorej przechowywany jest scoring danego poczatku pliku
#score of the begining of the peak
    Score_array = []
#Tablica w ktorej przechowywana jest informacja o numerze pliku z ktorego pochodzi numer
#(poczatka lub konca peaku)    
#file number of the begining peaks
    Begin_End_file_numbers_array = []
    
#Tablica w ktorej sprawdzamy czy dana wartosc to poczatek czy koniec kolumny
#check if a value is start or stop
    Start_Stop_array = []
    
    start = "start"
    stop = "stop"
    
#Liczba plikow, ktore udalo sie odczytac
#number of read files
    Read_files_number = len(line)
    

    del unrecognized_file_break
    del unread_files_numbers
    del substring

    

#W tej petli dziele linie odczytane z plikow (pierwsze linie z danymi)
#i dane wpisuje do odpowiednich tabeli (zadeklarowanych powyzej)    
# In this work, loop lines are read from the file (first line of data)
# And the data put in the appropriate table (declared above)
    for i in range(0, Read_files_number):
        splited_results = line[i].split('\t')
        Begin_End_array.append(splited_results[0])
        Score_array.append(str(int(splited_results[1])))
        stop_score = str(int(splited_results[1]) * (-1))
        Begin_End_file_numbers_array.append(i)
        Start_Stop_array.append(start)
        stop = int(splited_results[0]) + step
        Begin_End_array.append(str(stop))
        Score_array.append(stop_score)
        Begin_End_file_numbers_array.append(i)
        Start_Stop_array.append(stop)
        
 #Tablica z ustawionymi markerami konca pliku na 1. Jak to zostanie przestawione na 0, 
#to plik nie bedzie juz dalej odczytywany.    
#end of file marker
    EOF_markers = []
    
    
    for i in range(0, Read_files_number):
        EOF_markers.append(1)

#Zmienna pokazujaca liczbe odczytywanych nadal plikow. Jak spadnie do 0, program sie zatrzyma.        
#number of files that still has data 
    End_of_Files_marker = Read_files_number
    k = 0

#Sortujemy pobrane za pierwszym razem wartosci    
    Begin_End_array, Begin_End_file_numbers_array, Score_array, Start_Stop_array = Sort_numbers(Begin_End_array, Begin_End_file_numbers_array, Score_array, Start_Stop_array)



#Ta petla wypelnia plik tymczasowy    
    while End_of_Files_marker > 0:
        

#Pobieram dane z kompletem informacji (start/stop i numer pliku) i usuwam ja od razu z tablicy        
        operation_number = Begin_End_array.pop(0)

        
#Zmienna przechowujaca numer pliku, ktorym sie aktualnie zajmujemy        
        operation_file_number = Begin_End_file_numbers_array.pop(0)

        operation_marker = str(int(Score_array.pop(0)))
#Zmienna przechowujaca informacje czy to jest poczatek czy koniec peaku
        operation_start_stop = Start_Stop_array.pop(0)        
        
#Jesli start, to po prostu wpisze do pliku podnoszac jednoczesnie liczbe zliczen        
        if operation_start_stop == "start":
            k = k + 1
            
#Jesli stop, to poza zminiejszeniem liczby zliczen (o jedna pozycje za wczesnie) musze jeszcze zaczytac
#kolejna dana z pliku, z ktorego wiersz wlasnie skonczylismy
        else:
            k = k - 1
            line[operation_file_number] = ""
            line[operation_file_number] = results_files_id[operation_file_number].readline()
            
            
#Sprawdzam czy to nie byl koniec pliku
            if line[operation_file_number] != "":
                splited_results = line[operation_file_number].split("\t")
                Begin_End_array.append(splited_results[0])
                Score_array.append(str(int(splited_results[1])))
                stop_score = str(int(splited_results[1]) * (-1))   
                Begin_End_file_numbers_array.append(operation_file_number)
                Start_Stop_array.append(start)
                stop = int(splited_results[0]) + step
                Begin_End_array.append(str(stop))
                Score_array.append(stop_score)
                Begin_End_file_numbers_array.append(operation_file_number)
                Start_Stop_array.append(stop)
                Begin_End_array, Begin_End_file_numbers_array, Score_array, Start_Stop_array = Sort_numbers(Begin_End_array, Begin_End_file_numbers_array, Score_array, Start_Stop_array)                
            else:
                EOF_markers[operation_file_number] = 0
                End_of_Files_marker = 0
                for i in range(0, Read_files_number):
                    End_of_Files_marker = End_of_Files_marker + EOF_markers[i]
#Wpisuje dane do pliku tymczasowego                     
        if operation_marker != "0":
            temp_wig_compare_id.write(operation_number + "\t" + operation_marker + "\t" + str(k) + "\n")
#        print(operation_number + "\t" + operation_marker + "\t" + str(k))

#Zamykam wszystkie pliki wejsciowe
    del line
    for i in range(0, Read_files_number):              
        results_files_id[i].close()
        
    del results_files_id
    del EOF_markers
    del Begin_End_array
    del Score_array
    del Begin_End_file_numbers_array
    del splited_results

#Zamykam plik tymczasowy
    temp_wig_compare_id.close()

    return temp_wig_compare_file


