import pandas as pd
import numpy as np
import os
import sys

def get_keyres(ligs, file_list, output_file, **kwargs):

    #file_list=glob.glob('*.csv')
    pre_kres_file = kwargs.get('pre_kres_file', None)

    if pre_kres_file is not None:
        if os.path.exists(pre_kres_file):
            df1 = pd.read_csv(pre_kres_file, index_col = 'molecule')
            header_list = df1.columns.tolist()
        else:
            print('ERROR: %s does not exist.'%pre_kres_file)
            sys.exit(1)

    df_data= pd.DataFrame()
    #df = pd.read_csv('file.txt', sep=' ', header=None)
    fname=[]

    for i in file_list:
        if os.path.exists(i):
            df = pd.read_csv(i, sep='\s+', header=None)[12]
            #fname.append(os.path.basename(i))
            df_data = pd.concat([df_data, df], axis = 1, ignore_index=True)
        else:
            print('ERROR: %s does not exist.'%i)
            raise SystemExit
        

    df_data = df_data.T

    fname_new = []

    for name in ligs['id']:
        #temp = name.split('.', 1)[0]
        fname_new.append(name)

    df_data['molecule'] = fname_new

    df_data = df_data.set_index('molecule')

    df_data.columns = df_data.columns + 1

    #df_data

    drop_list = []

    if pre_kres_file is None:
        for res in df_data.columns:
            if min(df_data[res]) > -0.1:
                drop_list.append(res)
    else:
        for res in df_data.columns:
            if str(res) not in header_list:
                drop_list.append(res)

    df_keyres=df_data.drop(drop_list, axis = 1)

    df_keyres.to_csv(output_file, header = True)
