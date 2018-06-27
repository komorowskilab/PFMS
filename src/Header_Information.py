'''
Created on May 21, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.0

'''


def Get_Header_Info(file):
    
    file_id = open(file, "r")
    line = file_id.readline()
    substring = line[0:1]
    substring2 = line[0:3]
    unrecognized_file_break = 0
    status = 0
    step = -1
    chromosome = ""
    
    while substring.isdigit() == False and substring2 != "chr" and substring2 != "Chr" and substring2 != "CHR":
        arguments = line.split(" ")
        for j in range(0, len(arguments)):
            values = arguments[j].split("=")
            if values[0] == "span" or values[0] == "step":
#Zmienna step jest bardzo wazna. Przechowuje wartosc kroku pliku wiggle VariableStep 
#variable setp contains the value of the step.                       
                step = int(values[1])
            elif values[0] == "chrom":
                chromosome = values[1]
                                   
        line = file_id.readline()
        substring = line[0:1]
        substring2 = line[0:3]
        unrecognized_file_break = unrecognized_file_break + 1
        if unrecognized_file_break == 19:
            status = 1
          
    file_id.close()
      
    return status, step, chromosome
