'''
Created on May 15, 2010

@author: kruczyk
@organization: LCB centre for BioInformatics at Uppsala University
@version: 1.0

'''
def Fixed_to_Variable(input_file, label, step):
    span = str(step)
    FixedStep_file_id = open(input_file, "r")
    VariableStep_file_path = "./" + label + "_realwiggle.wig"
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
    

    
    score = int(float(list[0]))
    stop = start + FixedStep
    prev_stop = stop
    temp_step = step
    temp_score = 0
    till_step = 0
    
    while line != "":
        
        if till_step == 0:
            VariableStep_file_id.write(str(start) + "\t")
        
        if start + temp_step < stop:
            temp_score = (temp_score * till_step + score * temp_step)/step
            VariableStep_file_id.write(str(temp_score) + "\n")
            start = start + temp_step
            temp_step = step
            temp_score = 0
            till_step = 0
            
        elif start + temp_step == stop:
            temp_score = (temp_score * till_step + score * temp_step)/step
            VariableStep_file_id.write(str(temp_score) + "\n")
            temp_step = step
            temp_score = 0
            till_step = 0
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
                                    for i in range(0, len(list)):
                                        arguments = list[i].split("=")
                                        if arguments[0] == "chrom":
                                            temp_chromosome = arguments[1]
                                        
                                else:
                                    End_of_File = True
                                    VariableStep_file_id.write(str(temp_score) + "\n")
                                    break
                                    
                                if End_of_File == False:
                                    for i in range(0, len(list)):
                                        arguments = list[i].split("=")
                                        if arguments[0] == "start":
                                            start = int(arguments[1])
                                        elif arguments[0] == "step":
                                            FixedStep = int(arguments[1])
                                    
                    
                    stop = start + FixedStep
                    line = FixedStep_file_id.readline()
                    list = line.split(' ')

                elif list[0][0 : 1].isdigit() == True:
                    start = start + temp_step
                    stop = start + FixedStep
                    
                score = int(float(list[0]))
                

                
        else:
            temp_score = (temp_score * till_step + score * temp_step)/(till_step + temp_step)
            till_step = till_step + stop - start + 1
            temp_step = start + temp_step - stop
            prev_stop = stop
            line = FixedStep_file_id.readline()
            
            if line != "":
                list = line.split(' ')
                End_of_File = False
                
                if list[0][0 : 1].isdigit() == False:
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
                                    for i in range(0, len(list)):
                                        arguments = list[i].split("=")
                                        if arguments[0] == "chrom":
                                            temp_chromosome = arguments[1]
                                        
                                else:
                                    End_of_File = True
                                    VariableStep_file_id.write(str(temp_score) + "\n")
                                    break
                                    
                                if End_of_File == False:
                                    for i in range(0, len(list)):
                                        arguments = list[i].split("=")
                                        if arguments[0] == "start":
                                            start = int(arguments[1])
                                        elif arguments[0] == "step":
                                            FixedStep = int(arguments[1])

                    
                    line = FixedStep_file_id.readline()
                    list = line.split(' ')

                elif list[0][0 : 1].isdigit() == True:
                    start = start + temp_step
                    
                if End_of_File == False:
                    stop = start + FixedStep
                    score = int(float(list[0]))

            
                if start - prev_stop > 1 and prev_stop + temp_step >= start:
                    temp_score = temp_score * till_step
                    till_step = till_step + start - prev_stop - 1
                    temp_score = temp_score/till_step
                    temp_step = temp_step - (start - prev_stop - 1)
                    prev_stop = stop
                elif start - prev_stop > 1 and prev_stop + temp_step < start:
                    temp_score = temp_score * till_step / step
                    VariableStep_file_id.write(str(temp_score) + "\n")
                    temp_step = step
                    temp_score = 0
                    till_step = 0
            else:
                VariableStep_file_id.write(str(temp_score) + "\n")

        
    del list
    del arguments
    
    VariableStep_file_id.close()
    
    return VariableStep_file_path 

                              
    