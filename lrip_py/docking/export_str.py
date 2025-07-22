from schrodinger import structure
from schrodinger.utils import fileutils
from schrodinger.structutils import sort
import argparse
import json
#import pandas as pd

def export_str(input_file, out_file, out_property, count_file, lig_list, lig_format):
    # Post-docking structures and results processing stage
    # This module recongizes same structures (different docking poses) by their title (second row in mol2 files)
    # Can only read title from input without matching
    # Can only export mol2 structure sorted by glide score

    dic_title_cont = {}
    lig_list = json.loads(lig_list)
    # Determine the format of input file
    file_format = fileutils.get_structure_file_format(input_file)
    gscore = open(out_file, 'w')
    with structure.StructureReader(input_file, index=2) as reader:
        unq_st = []
        for st in reader:
            if st.title not in dic_title_cont:
                dic_title_cont[st.title] = 1
                unq_st.append(str(st.title))
                glide_score = st.property[out_property]
                gscore.write('%s_%s,%.4f\n'%(st.title, dic_title_cont[st.title], glide_score))
                with structure.StructureWriter('%s_%s.%s' % (st.title, dic_title_cont[st.title], lig_format)) as writer:
                    writer.append(st)
            else:
                dic_title_cont[st.title] += 1
                glide_score = st.property[out_property]
                gscore.write('%s_%s,%.4f\n'%(st.title, dic_title_cont[st.title], glide_score))
                with structure.StructureWriter('%s_%s.%s' % (st.title, dic_title_cont[st.title], lig_format)) as writer:
                    writer.append(st)

    with open(count_file, 'a') as cf:
        for nm in dic_title_cont:
            cf.write('%s,%d\n'%(nm,dic_title_cont[nm]))
    miss_list = []
    pos_flag = False
    for lig in lig_list:
        if lig not in unq_st:
            pos_flag = True
            miss_list.append(lig)
    if pos_flag:
        print('\nWARNING: Not all ligands have output poses. Missing ligands: ')
        for k in miss_list:
            print('%s '%k)
    
    with open('dock_fail_ligands.txt', 'w') as miss_file:
        if len(miss_list) == 0:
            miss_file.write('\n# No missing ligands found.')
        else:
            for lig in miss_list:
                miss_file.write("%s\n"%lig)


def sort_str(input_file, sort_key, out_file):
    sort.sort_file(input_file, sort_key, out_file_name=out_file)

def read_lig(lig_list_f):
    with open(lig_list_f, 'r') as f:
        lig_list = []
        for line in f:
            if len(line.strip()) > 0:
                lig_list.append(line.strip())
    
    return lig_list

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract Glide docking poses')
    parser.add_argument('-i', '--input', required=True, help='Post-docking _pv.maegz file')
    parser.add_argument('-o', '--out_file', required=True, help='Output .csv file')
    parser.add_argument('-p', '--out_property', required=True, help='The sequence sorted by the property, find the property key words on maestro website: https://learn.schrodinger.com/public/python_api/2023-2/api/schrodinger.application.glide.sort_utils.html?highlight=r_i_glide_gscore.')
    parser.add_argument('-c', '--count_file', required=True, help='File to count the number of output poses')
    parser.add_argument('-l', '--lig_list', required=False)
    parser.add_argument('-ll', '--lig_list_f', required=False)
    parser.add_argument('-f', '--lig_format', required=True)
    args = parser.parse_args()
    return args
if __name__ == '__main__':
    # input_file = "./test_docking_pv.maegz"
    # find properties from: https://learn.schrodinger.com/public/python_api/2023-2/api/schrodinger.application.glide.sort_utils.html?highlight=r_i_glide_gscore
    # sort_key = [('r_i_glide_gscore', sort.ASCENDING)]
    # out_property = 'r_i_glide_gscore'
    args            = parse_arguments()
    input_file      = args.input
    out_file        = args.out_file
    out_property    = args.out_property
    count_file      = args.count_file
    lig_list        = args.lig_list
    lig_list_f      = args.lig_list_f
    lig_format      = args.lig_format

    if lig_list:
        lig_list = lig_list
    elif lig_list_f:
        #lig_list_tmp = read_lig(lig_list_f)
        print("Reading ligand list from: %s"%lig_list_f)
        with open(lig_list_f, 'r') as f:
            lig_list_tmp = json.loads(f)
            lig_list = json.dumps(lig_list_tmp)

    export_str(input_file, out_file, out_property, count_file, lig_list, lig_format)
        

