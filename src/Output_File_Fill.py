'''
@attention: Writes the comparison results to the final output file
Created on May 27, 2010
@author: kruczyk
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 1 Mar 2011, , some fixes
@version: 1.0
'''

def Fill_Output_File(temp_wig_compare_file, label, chromosome, step):
    
    output_file_name = "./" + label+"_Results.wig"
    output_file_id = open(output_file_name, "w")
    
    header = "track type=wiggle_0 name=\"" + label + "_full_PFMeataserver\" description=\"Peak Finder Metaserver (All Peaks) results\" color=50,50,150 visibility=2\nvariableStep chrom=" + chromosome + " span=" + str(step) + "\n"
    output_file_id.write(header)
    
    #Otwieram plik tymczasowy do czytania   
    
    temp_wig_compare_id = open(temp_wig_compare_file, "r")

#Odczytuje pierwsza linie z pliku tymczasowego    
    line = temp_wig_compare_id.readline()


    
# W tej petli przeprowadzone jest glosowanie i uzupelniany jest plik wynikowy WIG
#the voting mechanism is implemented here, working on the temp file that was created by TEmp_fill_file and write it to wiggle file
    start_of_pile = 0
    stop_of_pile = 0
    actual_score = 0
    till_points = 0
    till_step = 0
    list = []

    while line != "":
        
        list = line.split('\t')
        read_start = int(list[0])
        read_score = int(list[1])

        
        if read_score == 0 and actual_score == 0:
            line = temp_wig_compare_id.readline()
            continue  
        

        if read_start <= stop_of_pile:
            till_points = till_points + (actual_score * (read_start - start_of_pile - till_step))
            till_step = read_start - start_of_pile
            actual_score = actual_score + read_score
        else:
            till_points = till_points + (actual_score * (step - till_step))
            rank = int(round(till_points / float(step)))
            if rank > 0:
                output_file_id.write(str(start_of_pile) + "\t")
                output_file_id.write(str(rank) + "\n")
                
            if actual_score > 0:
                start_of_pile = stop_of_pile + 1
                stop_of_pile = stop_of_pile + step
                while read_start > stop_of_pile:
                    output_file_id.write(str(start_of_pile) + "\t")
                    output_file_id.write(str(actual_score) + "\n")
                    start_of_pile = start_of_pile + step
                    stop_of_pile = stop_of_pile + step
                        
                till_step = read_start - start_of_pile
                till_points = actual_score * till_step
                actual_score = actual_score + read_score
                    
            else:
                start_of_pile = read_start
                stop_of_pile = read_start + step - 1
                actual_score = read_score
                till_step = 0
                till_points = 0
                        
            
        line = temp_wig_compare_id.readline()
        
        
    del list
    output_file_id.close()
    temp_wig_compare_id.close()
    
    
    return output_file_name    
            
            
                
        
        
        
        
        
        
