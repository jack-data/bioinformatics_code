import sys #Needed for accessing cmd line args
import re #Needed for regex to match various patterns
import os, glob #Needed for opening directories and parsing paths

#Script for reading FASTA files and calculating metrics on them.
#Can take individual FASTA files as arguments or take a directory containing FASTA files.
#output is a single text file, its name is an argument by user using -o flag.
#
#Text file is of format:
#
#Filename: <filename>
#Total Contig Length: <Total contig length>
#Total GC content: <Total GC content>
#GC percentage: <GC percentage>
#N50: <n50>
#N90: <n90>
#L50: <l50>
#       Header: <header>
#           Length: <contig length>
#           GC content: <GC content>
#       Header: <header>
#           Length: <contig length>
#           GC content: <GC content>
#       ...

#REGEX patterns
FLAG_PATTERN = '-[a-z$]' #Pattern for handling various flags at cmd line such as -h
FASTA_FILE_PATTERN = '[a-zA-Z0-9]+\.fa$' #Pattern for verifying that a file is of *.fa format (and excludes for example *.fa.exe format)
ANY_FILE_PATTERN = '\.+[a-zA-Z0-9]+'

#Message that is printed when the script is ran with -h flag (help menu)
HELP_MESSAGE = ('\nHELP MENU\nScript takes format of: <scriptname> <files/directories> <flags>\n'
'If a directory name contains a space it must be surrounded by quotation marks like so: \'dir ec tory\'\n\n'
'Flags include:\n\n'
'-h: Prints this menu (and exits script).\n\n'
'-o: (REQUIRED) Output to text file; requires an argument to follow which will be used as the filename. \n'
'\tFilename cannot be the same as an existing file.\n')

#Functions
def open_FASTA(file):
    #Open a FASTA file.
    #Returns as a dictionary of header:contig pairs.
    with open(file, "r") as file_handler: #File is read only, use with keyword to ensure file is closed afterwards.
        map = {} #Dictionary to store header:contigs
        header = ''
        contig = ''
        for line_number, line in enumerate(file_handler): #Use enumerate so that we don't have to load the entire file into memory.
            line = line.strip() #Remove newline characters from line
            if line[0] == '>' and header == '': #First header
                header = line
            else:
                if((line[0] != '>')): #Line is not a header
                    contig += line #Build contig
                else: #Line is a header
                    map[header] = contig #Store prev. header and contig as key:value
                    header = line #New header, continue
                    contig = '' #Reset contig
        map[header] = contig #Reached end of file so include the last key:value we had (as there are no more lines which are headers)
    return map #Return dictionary of headers:reads

def get_contig_length(contig):
    #Get the length of a contig.
    #Returns as an integer.
    return len(contig)

def get_GC_content(contig):
    #Get the GC content of a contig.
    #Returns as an integer.
    return contig.count('G') + contig.count('C') 

def get_total_GC_content(opened_FASTA):
    #Get the Total GC content of a FASTA file.
    #Use open_FASTA() to obtain opened_FASTA parameter.
    #Returns as an integer.
    contigs = list(opened_FASTA.values())
    return sum(get_GC_content(contig) for contig in contigs)

def get_total_contig_length(opened_FASTA):
    #Get the Total length of all contigs in a FASTA file.
    #Use open_FASTA() to obtain opened_FASTA parameter.
    #Returns as an integer. 
    contigs = list(opened_FASTA.values())
    return sum(len(contig) for contig in contigs)

def get_GC_percentage(GC, length):
    return (GC/length) * 100

def get_N_or_F_value(opened_FASTA, val, flag):
    #Get either the Nval or Fval of a FASTA file. 
    #For example to get the N50 value, run with val as 50 and 'N' as flag.
    #Use open_FASTA() to obtain opened_FASTA parameter.
    #Returns as an integer.
    contigs = list(opened_FASTA.values()) #Store just the contig values from opened file into a list
    contigs.sort(key=len, reverse=True) #Contigs list is now sorted in descending order of length

    #Iterate through to get combined total length of all contigs.
    total_length = get_total_contig_length(opened_FASTA)

    #Iterate through again, checking against total length to find appropriate metric.
    running_total = 0
    for contig in contigs:
        running_total += len(contig)
        if(running_total >= (total_length*(val/100))):
            if(flag == 'N'):
                return len(contig)
            elif(flag == 'F'):
                return contigs.index(contig)+1 #Have to add one to return value, as L50 value is 1-indexed rather than 0-indexed like lists are.
            else:
                raise Exception('Incorrect flag argument, must be \'N\' or \'F\'')


def handle_arguments():
    #Handles arguments given at command line.
    #Reads each argument and for each argument checks if it is a flag, file, directory or output filename.
    output = {}
    output_flag = False

    if(len(sys.argv) > 1): #Check that arguments have been given
        if('-h' in sys.argv): #If -h is found as any argument will only print help message. (saves wasted computation if for example: "<script> <large directory> -h")
            print(HELP_MESSAGE)
            return

        for arg in sys.argv:

            if(arg == sys.argv[0]): #Ignore script name as argument
                continue
            elif(output_flag):
                if(re.search(ANY_FILE_PATTERN,arg)): #If no file format is given, '.txt' will be given by default
                    output_file_name = arg
                else:
                    output_file_name = arg + '.txt'

                with open(output_file_name, 'x') as output_file: #Create output file
                    for opened_FASTA in output:
                        total_contig_length = get_total_contig_length(output[opened_FASTA])
                        total_GC_content = get_total_GC_content(output[opened_FASTA])
                        GC_percentage = get_GC_percentage(total_GC_content, total_contig_length)
                        N50_value = get_N_or_F_value(output[opened_FASTA], 50, 'N')
                        N90_value = get_N_or_F_value(output[opened_FASTA], 90, 'N')
                        L50_value = get_N_or_F_value(output[opened_FASTA], 50, 'F')

                        output_file.write('File: '+opened_FASTA+'\n')
                        output_file.write('Total Contig Length: '+str(total_contig_length)+'\n')
                        output_file.write('Total GC Content: '+str(total_GC_content)+'\n')
                        output_file.write('GC percentage: '+str(GC_percentage)+'%\n')
                        output_file.write('N50: '+str(N50_value)+'\n')
                        output_file.write('N90: '+str(N90_value)+'\n')
                        output_file.write('L50: '+str(L50_value)+'\n')

                        for header in output[opened_FASTA]:
                            output_file.write('\tHeader: '+header+'\n')
                            output_file.write('\t\tLength: '+str(get_contig_length(output[opened_FASTA][header]))+'\n')
                            output_file.write('\t\tGC Content: '+str(get_GC_content(output[opened_FASTA][header]))+'\n')
                output_flag = False
                return
            elif(re.search(FLAG_PATTERN, arg)): #Is argument a flag that needs to look at next argument? e.g. -o
                match(arg):
                    case '-o':
                        #Ensures that next argument will be read as the output text file
                        output_flag = True
            elif(re.search(FASTA_FILE_PATTERN, arg)): #Argument not flag, is it a fasta file?
                output[arg] = open_FASTA(arg)
            elif(glob.glob(os.path.join(arg, '*.fa'))):
                #Returns true if argument is a path to a folder containing fasta files
                for file in glob.glob(os.path.join(arg, '*.fa')):
                    output[file] = open_FASTA(file)
            else:
                raise Exception('\nArguments appear invalid. Try running with -h flag for help.')
    else:
        raise Exception('\nArguments must be given. Try running with -h flag for help.')

if __name__ == '__main__':
    handle_arguments()