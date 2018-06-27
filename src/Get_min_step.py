'''
Created on May 28, 2010
@author: kruczyk
'''

def Get_step(file):
    
    file_id = open(file, "r")
    line = file_id.readline()
    substring = line[0:1]
    substring2 = line[0:3]
    unrecognized_file_break = 0
    status = 0
    step = -1
    
    while substring.isdigit() == False and substring2 != "chr" and substring2 != "Chr" and substring2 != "CHR":
        arguments = line.split(" ")
        for j in range(0, len(arguments)):
            values = arguments[j].split("=")
            if values[0] == "span" or values[0] == "step":
#Zmienna step jest bardzo wazna. Przechowuje wartosc kroku pliku wiggle VariableStep                        
                step = int(values[1])

                                   
        line = file_id.readline()
        substring = line[0:1]
        substring2 = line[0:3]
        unrecognized_file_break = unrecognized_file_break + 1
        if unrecognized_file_break == 19:
            status = 1
          
    file_id.close()
      
    return status, step

def Get_min_step(files):
    
    min_step = -1
    
    for i in range(0, len(files)):
        status, temp_step = Get_step(files[i])

        if status == 0:

                
            if temp_step > 0:
                if min_step <= 0:
                    min_step = temp_step
                else:
                    if temp_step < min_step:
                        min_step = temp_step 
                        
    return min_step        
    