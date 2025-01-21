import os
import sys

def extract_ie(input_file, comp_info, output_file):
    if os.path.exists(input_file):
        pass
    else:
        print('ERROR: %s doesn\'t exist, exit.'%input_file)
        sys.exit(1)
    output = open(output_file, 'w')
    with open(input_file) as input:
        for line in input:
            split_line = line.split()
            if split_line[0] != '#':
                if int(split_line[0]) >= int(comp_info['ligsrt']) and int(split_line[0]) <= int(comp_info['ligend']) and int(split_line[0]) != int(split_line[1]):
                    output.write(line)
