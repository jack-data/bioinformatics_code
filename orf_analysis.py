import sys
import re

ORF_pattern = '[0-9]+-[0-9]+'

def sortORFlengths(file):
    orf_lengths = {}

    with open(file, "r") as file_handler:
        for line_number, line in enumerate(file_handler): #Use enumerate so that we don't have to load the entire file into memory.
            if(line[0] == '>'): #Only operate on lines with headers
                header = line.strip() #Remove newlines

                orf_indexes = re.search(ORF_pattern, header).span() #Gets the index of for example '148-3535' in the header

                orf_start_and_end = line[orf_indexes[0]:orf_indexes[1]].split('-') #Split into list of first number (beginning pos. of ORF) and last (end pos.)
                orf_start = orf_start_and_end[0]
                orf_end = orf_start_and_end[1]
                        
                #Minus the substring from first part of string to get ORF length
                orf_length = int(orf_end) - int(orf_start)

                orf_lengths[header] = orf_length #Store in dictionary

    orf_lengths = sorted(orf_lengths.items(), key=lambda x: x[1], reverse=True) #Sort dictionary into list of tuples in descending order by size of orf_length

    with open((file+' sorted orfs.txt'), 'x') as output_file: #Write the output to a file
        for orf_length in orf_lengths:
            output_file.write(str(orf_length[0])+' Length: '+str(orf_length[1])+'\n')

if __name__ == '__main__':
    sortORFlengths(sys.argv[1])