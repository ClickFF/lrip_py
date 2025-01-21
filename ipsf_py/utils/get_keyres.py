import pandas as pd
import numpy as np
import os

def get_keyres(ligs, file_list, output_file):

    #file_list=glob.glob('*.csv')

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

    for res in df_data.columns:
        if min(df_data[res]) > -0.1:
            drop_list.append(res)

    df_keyres=df_data.drop(drop_list, axis = 1)

    #df_keyres

    df_keyres.to_csv(output_file, header = True)
