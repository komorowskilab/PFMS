'''
Created on May 16, 2010

@author: kruczyk
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 14 Jan 2011, some fixes
@version: 1.0

'''



def Variable_to_Variable(input_file, label, step):
    input_VariableStep_file_id = open(input_file, "r")
    line = input_VariableStep_file_id.readline()
    list = line.split(' ')
    lines = []
    File_is_ok = False
    
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
                    
        

        if original_step == step:
            File_is_ok = True
            break
        else:
            lines.append(line)
            line = input_VariableStep_file_id.readline()
            list = line.split(' ')

    
    if File_is_ok == True:
        del list
        del lines
        del arguments
        print "File is OK!"
        output_VariableStep_file = input_file
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
        score = int(float(list[1]))
        prev_stop = stop
        temp_step = step
        temp_score = 0
        till_step = 0
    
    
        while line != "":
            if till_step == 0:
                output_VariableStep_file_id.write(str(start) + "\t")

        
            if start + temp_step - 1 < stop:
                temp_score = (temp_score * till_step + score * temp_step)/step
                output_VariableStep_file_id.write(str(temp_score) + "\n")
                start = start + temp_step
                temp_step = step
                temp_score = 0
                till_step = 0

        
            elif start + temp_step - 1 == stop:
                temp_score = (temp_score * till_step + score * temp_step)/step
                output_VariableStep_file_id.write(str(temp_score) + "\n")
                temp_step = step
                temp_score = 0
                till_step = 0
                line = input_VariableStep_file_id.readline()
            
                if line != "":
                    list = line.split('\t')
                    start = int(list[0])
                    stop = start + original_step - 1
                    score = int(float(list[1]))
                  
                
            else:
                temp_score = (temp_score * till_step + score * temp_step)/(till_step + temp_step)
                till_step = till_step + stop - start + 1
                temp_step = start + temp_step - stop
                prev_stop = stop
                line = input_VariableStep_file_id.readline()
            
                if line != "":
                    list = line.split('\t')
                    start = int(list[0])
                    stop = start + original_step - 1
                    score = int(float(list[1]))
                    
                    if start - prev_stop > 1 and prev_stop + temp_step >= start:
                        temp_score = temp_score * till_step
                        till_step = till_step + start - prev_stop - 1
                        temp_score = temp_score/till_step
                        temp_step = temp_step - (start - prev_stop - 1)
                        prev_stop = stop
                    elif start - prev_stop > 1 and prev_stop + temp_step < start:
                        temp_score = temp_score * till_step / step
                        output_VariableStep_file_id.write(str(temp_score) + "\n")
                        temp_step = step
                        temp_score = 0
                        till_step = 0
                else:
                    output_VariableStep_file_id.write(str(temp_score) + "\n")
                    
        del list
        output_VariableStep_file_id.close()

    
    return output_VariableStep_file
            
            