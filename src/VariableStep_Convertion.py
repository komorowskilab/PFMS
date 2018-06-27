'''
Created on May 25, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.0

'''

import shutil

def BedGraph_to_Variable(input_file, label, step):
    span = str(step)
    BedGraph_file_id = open(input_file, "r")
    Wiggle_file_path = "./" + label + "_realwiggle_bedgraph.wig"
    Wiggle_file_id = open(Wiggle_file_path, "w")
    Wiggle_file_id.write("track type=wiggle_0 name=" + label + " ")
    line = BedGraph_file_id.readline()
    list = line.split(' ')
    while list[0][0:3] != "chr" and list[0][0:3] != "Chr" and list[0][0:3] != "CHR":
        for i in range(0, len(list)):
            arguments = list[i].split("=")
            
            if arguments[len(arguments) - 1][len(arguments[len(arguments) - 1]) - 1 : len(arguments[len(arguments) - 1])] == "\n":
                list[i] = list[i][0 : len(list[i]) - 1]            
            
            if arguments[0] in ("description", "visibility", "color", "itemRgb", "colorByStrand", "useScore", "group", "priority", "db", "offset", "url", "htmlUrl"):
                Wiggle_file_id.write(list[i] + " ")                      
        line = BedGraph_file_id.readline()
        list = line.split(' ')
    Wiggle_file_id.write("\n")
    Wiggle_file_id.write("variableStep chrom=" + list[0] + " span=" + span + "\n")
    
    start = int(list[1])
    stop = int(list[2]) - 1
    score = int(round(float(list[3])))
    prev_stop = stop
    temp_step = step
    temp_score = 0
    till_step = 0
    
    while line != "":
        if till_step == 0:
            Wiggle_file_id.write(str(start) + "\t")
        
        if start + temp_step - 1 < stop:
            temp_score = (temp_score + score * temp_step)
            out_score = int(round(temp_score/float(step)))
            Wiggle_file_id.write(str(out_score) + "\n")
            start = start + temp_step
            temp_step = step
            temp_score = 0
            till_step = 0

        
        elif start + temp_step - 1 == stop:
            temp_score = (temp_score + score * temp_step)
            out_score = int(round(temp_score/float(step)))
            Wiggle_file_id.write(str(out_score) + "\n")
            temp_step = step
            temp_score = 0
            till_step = 0
            line = BedGraph_file_id.readline()
            
            if line != "":
                list = line.split(' ')
                start = int(list[1])
                stop = int(list[2]) - 1
                score = int(round(float(list[3])))
                
        else:
            old_till_step = till_step
            till_step = till_step + stop - start + 1
            temp_step = start + temp_step - stop - 1
            temp_score = (temp_score + score * (till_step - old_till_step))
            prev_stop = stop
            line = BedGraph_file_id.readline()
            
            if line != "":
                list = line.split(' ')
                start = int(list[1])
                stop = int(list[2]) - 1
                score = int(round(float(list[3])))

            
                if start - prev_stop > 1 and prev_stop + temp_step >= start:
                    
                    till_step = till_step + start - prev_stop - 1
                    temp_step = temp_step - (start - prev_stop - 1)
                    prev_stop = stop
                    
                elif start - prev_stop > 1 and prev_stop + temp_step < start:
                    
                    out_score = int(round(temp_score/float(step)))
                    Wiggle_file_id.write(str(out_score) + "\n")
                    temp_step = step
                    temp_score = 0
                    till_step = 0
            else:
                out_score = int(round(temp_score/float(step)))
                Wiggle_file_id.write(str(out_score) + "\n")
    
    del list
    Wiggle_file_id.close()
    
    return Wiggle_file_path    
            

def Fixed_to_Variable(input_file, label, step):
    span = str(step)
    FixedStep_file_id = open(input_file, "r")
    VariableStep_file_path = "./" + label + "_realwiggle_fixed.wig"
    VariableStep_file_id = open(VariableStep_file_path, "w")
    VariableStep_file_id.write("track type=wiggle_0 name=" + label + " ")
    line = FixedStep_file_id.readline()
    list = line.split(' ')
    
    while list[0][0 : 1].isdigit() == False:
#        print line[0 : len(list[0]) - 1]
#        print line[0 : len(list[0]) - 1].isdigit()
#        print len(line[0 : len(list[0]) - 1])
#        print ~line[0 : len(list[0]) - 1].isdigit()
#        print list[0].isdigit()
#        print (~list[0].isdigit() and ~line[0 : len(list[0]) - 1].isdigit())
        for i in range(0, len(list)):
            arguments = list[i].split("=")
            
            if arguments[len(arguments) - 1][len(arguments[len(arguments) - 1]) - 1 : len(arguments[len(arguments) - 1])] == "\n":
                list[i] = list[i][0 : len(list[i]) - 1]
            
            
            if arguments[0] in ("description", "visibility", "color", "itemRgb", "colorByStrand", "useScore", "group", "priority", "db", "offset", "url", "htmlUrl"):
                VariableStep_file_id.write(list[i] + " ")
            elif arguments[0] == "chrom":
                chromosome = arguments[1]
            elif arguments[0] == "start":
                start = int(arguments[1]) 
            elif arguments[0] == "step":
                FixedStep = int(arguments[1])
                                      
        line = FixedStep_file_id.readline()
        list = line.split(' ')

    VariableStep_file_id.write("\n")
    VariableStep_file_id.write("variableStep chrom=" + chromosome + " span=" + span + "\n")
    

    
    score = int(round(float(list[0])))
    stop = start + FixedStep - 1
    prev_stop = stop
    temp_step = step
    temp_score = 0
    till_step = 0
    
    while line != "":
        

        if till_step == 0:
            VariableStep_file_id.write(str(start) + "\t")

        
        if start + temp_step - 1 < stop:
            temp_score = (temp_score + score * temp_step)
            out_score = int(round(temp_score/float(step)))
            VariableStep_file_id.write(str(out_score) + "\n")
            start = start + temp_step
            temp_step = step
            temp_score = 0
            till_step = 0
            
        elif start + temp_step - 1 == stop:
            
            temp_score = (temp_score + score * temp_step)
            out_score = int(round(temp_score/float(step)))
            VariableStep_file_id.write(str(out_score) + "\n")
            line = FixedStep_file_id.readline()
            
            if line != "":
                list = line.split(' ')
                End_of_File = False
                if list[0][0 : 1].isdigit() == False:
                    old_start = start
                    old_FixedStep = FixedStep

                    for i in range(0, len(list)):
                        arguments = list[i].split("=")
                        if arguments[0] == "start":
                            start = int(arguments[1])
                        elif arguments[0] == "step":
                            FixedStep = int(arguments[1])
                        elif arguments[0] == "chrom" and arguments[1] != chromosome:
                            start = old_start
                            FixedStep = old_FixedStep
                            temp_chromosome = arguments[1]                            
                            while temp_chromosome != chromosome:
                                line = FixedStep_file_id.readline()
                                
                                if line != "":
                                    list = line.split(' ')
                                    for j in range(0, len(list)):
                                        arguments = list[j].split("=")
                                        if arguments[0] == "chrom":
                                            temp_chromosome = arguments[1]
                                        
                                else:
                                    End_of_File = True
                                    break
                                    
                                if End_of_File == False:
                                    for j in range(0, len(list)):
                                        arguments = list[j].split("=")
                                        if arguments[0] == "start":
                                            start = int(arguments[1])
                                        elif arguments[0] == "step":
                                            FixedStep = int(arguments[1])
                            
                            if End_of_File == True:
                                break        
                    
                    stop = start + FixedStep - 1
                    line = FixedStep_file_id.readline()
                    list = line.split(' ')

                elif list[0][0 : 1].isdigit() == True:
                    start = start + temp_step
                    stop = start + FixedStep - 1
                    
                if End_of_File == False:
                    score = int(round(float(list[0])))
            
            temp_step = step
            temp_score = 0
            till_step = 0    

                
        else:
            old_till_step = till_step
            till_step = till_step + stop - start + 1
            temp_step = start + temp_step - stop - 1
            temp_score = (temp_score + score * (till_step - old_till_step))
            prev_stop = stop
            line = FixedStep_file_id.readline()
            
            if line != "":
                list = line.split(' ')
                End_of_File = False
                
                if list[0][0 : 1].isdigit() == False:
                    old_start = start
                    old_FixedStep = FixedStep
                    for i in range(0, len(list)):
                        if End_of_File == True:
                            break
                        
                        arguments = list[i].split("=")
                        if arguments[0] == "start":
                            start = int(arguments[1])
                        elif arguments[0] == "step":
                            FixedStep = int(arguments[1])
                        elif arguments[0] == "chrom" and arguments[1] != chromosome:
                            start = old_start
                            FixedStep = old_FixedStep
                            temp_chromosome = arguments[1]                            
                            while temp_chromosome != chromosome:
                                line = FixedStep_file_id.readline()
                                
                                if line != "":
                                    list = line.split(' ')
                                    for j in range(0, len(list)):
                                        arguments = list[j].split("=")
                                        if arguments[0] == "chrom":
                                            temp_chromosome = arguments[1]
                                        
                                else:
                                    End_of_File = True
                                    out_score = int(round(temp_score/float(step)))
                                    VariableStep_file_id.write(str(out_score) + "\n")
                                    break
                                    
                                if End_of_File == False:
                                    for j in range(0, len(list)):
                                        arguments = list[j].split("=")
                                        if arguments[0] == "start":
                                            start = int(arguments[1])
                                        elif arguments[0] == "step":
                                            FixedStep = int(arguments[1])
                            
                            
                            if End_of_File == True:
                                break                       
                    
                    line = FixedStep_file_id.readline()
                    list = line.split(' ')

                elif list[0][0 : 1].isdigit() == True:
                    start = stop + 1
                    
                if End_of_File == False:
                    stop = start + FixedStep - 1
                    score = int(round(float(list[0])))

            
                if start - prev_stop > 1 and prev_stop + temp_step >= start:
                    till_step = till_step + start - prev_stop - 1
                    temp_step = temp_step - (start - prev_stop - 1)
                    prev_stop = stop
                    
                elif start - prev_stop > 1 and prev_stop + temp_step < start:
                    
                    out_score = int(round(temp_score/float(step)))
                    VariableStep_file_id.write(str(out_score) + "\n")
                    temp_step = step
                    temp_score = 0
                    till_step = 0
            else:
                out_score = int(round(temp_score/float(step)))
                VariableStep_file_id.write(str(out_score) + "\n")

        
    del list
    del arguments
    
    VariableStep_file_id.close()
    
    return VariableStep_file_path 


def Variable_to_Variable(input_file, label, step):
    input_VariableStep_file_id = open(input_file, "r")
    line = input_VariableStep_file_id.readline()
    list = line.split(' ')
    lines = []
    File_is_ok = False
    
#Ustawiam sie na pierwszej linii za naglowkiem    
#Zczytuje tez nazwe pliku i krok
    
    
    while list[0][0 : 1].isdigit() == False:
        line = ""
        original_step = ""
        for i in range(0, len(list)):
            arguments = list[i].split("=")

            if arguments[0] == "span":
                original_step = int(arguments[1])
                if list[i][len(list[i]) - 1 : len(list[i])] == "\n":
                    list[i] = "span=" + str(step) + "\n"
                else:
                    list[i] = "span=" + str(step)
                    
            elif arguments[0] == "name":
                if list[i][len(list[i]) - 1 : len(list[i])] == "\n":                
                    list[i] = "name=" + label + "\n"
                else:
                    list[i] = "name=" + label
            
            if line == "":
                line = list[i]
            else:
                line = line + " " + list[i]
                    
        
#Jesli krok jest taki sam, na jaki chcemy zmienic, to wychodzimy
        if original_step == step:
            File_is_ok = True
            break
#Jesli jest inaczej, to linie naglowka wpisuje do lines i nowa linie wczytuje do zmiennej line
        else:
            lines.append(line)
            line = input_VariableStep_file_id.readline()
            list = line.split(' ')

#Jesli pliku nie trzeba przeksztalcac, to sprzatamy po sobie i wychodzimy    
    if File_is_ok == True:
        del list
        del lines
        del arguments
        print "File is OK!"
        shutil.copy(input_file, "ranked_"+input_file)
        output_VariableStep_file = "ranked_"+input_file
#W przeciwnym razie otwieramy plik wynikowy i zaczynamy zabawe
    else:
        output_VariableStep_file = "./" + label + "_new_step.wig"
        output_VariableStep_file_id = open(output_VariableStep_file, "w")
        for i in range(0, len(lines)):
            temp_word = lines.pop(0)
            output_VariableStep_file_id.write(temp_word)
       
        del lines

        list = line.split('\t') 
        start = int(list[0])
        stop = start + original_step - 1
        score = int(round(float(list[1])))
        prev_stop = stop
        temp_step = step
        temp_score = 0
        till_step = 0
        
    
    
        while line != "":
            
            if till_step == 0:
                output_VariableStep_file_id.write(str(start) + "\t")
            

#Jesli nasz krok miesci sie w obrebie poprzedniego odczytu wchodzimy tutaj        
            if start + temp_step - 1 < stop:
                temp_score = (temp_score + score * temp_step)
                out_score = int(round(temp_score/float(step)))
                output_VariableStep_file_id.write(str(out_score) + "\n")
                start = start + temp_step
                temp_step = step
                temp_score = 0
                till_step = 0

#Jesli nasz krok konczy sie w tym miejscy, co poprzedni odczyt to wchodzimy tutaj
            elif start + temp_step - 1 == stop:
#Obliczamy temp_score - czyli iloczyn dotychczasowych par bazowych i ich scoringow 
                temp_score = (temp_score + score * temp_step)
                out_score = int(round(temp_score/float(step)))                  
                output_VariableStep_file_id.write(str(out_score) + "\n")
                temp_step = step
                temp_score = 0
                till_step = 0
                line = input_VariableStep_file_id.readline()
            
                if line != "":
                    list = line.split('\t')
                    start = int(list[0])
                    stop = start + original_step - 1
                    score = int(round(float(list[1])))
                  
                
            else:
                old_till_step = till_step
                till_step = till_step + stop - start + 1
                temp_step = start + temp_step - stop - 1
                temp_score = (temp_score + score * (till_step - old_till_step))
                prev_stop = stop
                line = input_VariableStep_file_id.readline()
            
                if line != "":
                    list = line.split('\t')
                    start = int(list[0])
                    stop = start + original_step - 1
                    score = int(round(float(list[1])))
                    
                    if start - prev_stop > 1 and prev_stop + temp_step >= start:
 
                        till_step = till_step + start - prev_stop - 1
                        temp_step = temp_step - (start - prev_stop - 1)
                        prev_stop = stop
                        
                    elif start - prev_stop > 1 and prev_stop + temp_step < start:

                        out_score = int(round(temp_score/float(step)))
                        output_VariableStep_file_id.write(str(out_score) + "\n")
                        temp_step = step
                        temp_score = 0
                        till_step = 0
                else:
                    out_score = int(round(temp_score/float(step)))
                    output_VariableStep_file_id.write(str(out_score) + "\n")
                    
        del list
        output_VariableStep_file_id.close()

    
    return output_VariableStep_file


def Convert_to_VariableStep(file, label, step):

    file_id = open(file, "r")
    line = file_id.readline()

    
    type = ""
    bedGraph = ["bedGraph", "BedGraph", "BEDGRAPH", "bedgraph"]
    variableStep = ["variableStep", "VariableStep", "VARIABLESTEP", "variablestep"]
    fixedStep = ["fixedStep", "FixedStep", "FIXEDSTEP", "fixedstep"]
    
    
    while type == "" and line != "" and line[0:1].isdigit() == False and line[0:3] != "chr" and line[0:3] != "Chr" and line[0:3] != "CHR" and line[0 : len(label)] != label:
        arguments = line.split(" ")
        for i in range(0, len(arguments)):
            argument = arguments[i].split("=")
            if argument[0] == "type":
                if argument[1] in bedGraph:
                    type = "bedGraph"
            elif argument[0] in variableStep:
                type = "variableStep"

            elif argument[0] in fixedStep:
                type = "fixedStep"
        line = file_id.readline()
    
    file_id.close()
    
    if type == "bedGraph":
        output_file = BedGraph_to_Variable(file, label, step)
        status = 0
        
    elif type == "variableStep":
        output_file = Variable_to_Variable(file, label, step)
        status = 0
    
    elif type == "fixedStep":
        output_file = Fixed_to_Variable(file, label, step)
        status = 0
    
    else:
        status = 1
        
    del bedGraph
    del fixedStep
    del variableStep
    
        
    return output_file, status
