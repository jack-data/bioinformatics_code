import sys
import re

#Script used for reshaping data from FASTA file.

HEADER_PATTERN = '>+[a-z]+_[0-9]+:+[0-9]+-[0-9]+\(+[+-]+\)'
#Example header:
#>scaffold_1693:3657-8985(-)

def get_n_lines(n, file):
    lines = []
    with open(file, "r") as file_handler:
        for line_number, line in enumerate(file_handler):
            if line_number == n:
                break

            line = line.strip()
            lines.append(line)
    return lines

def find_n_headers_in_FASTA(n, orfs_file, file_to_find_headers):
    n_lines = get_n_lines(n, orfs_file)
    headers = []
    for line in n_lines:
        indexes = re.search(HEADER_PATTERN, line).span()
        headers.append(line[indexes[0]:indexes[1]])

    found_headers_with_contigs = {}

    for header in headers:
        with open(file_to_find_headers, "r") as file_handler:
            contig = ''
            header_found = False
            
            for line_number, line in enumerate(file_handler):
                line = line.strip()
                if(header in line):
                    header_found = True
                    continue
                if(header_found):
                    if(line[0] != '>'):
                        contig += line
                    elif(line[0] == '>'):
                        found_headers_with_contigs[header] = contig
                        contig = ''
                        header_found = False

    return found_headers_with_contigs


if __name__ == '__main__':
    found = find_n_headers_in_FASTA(int(sys.argv[1]), sys.argv[2], sys.argv[3]) #Map of header:contigs
    found = list(found.items()) #Convert to list

    output_file_name = sys.argv[4]

    with open(output_file_name, "x") as output_file:
            for item in found:
                output_file.write(str(item[0])+'\n')
                output_file.write(str(item[1])+'\n')