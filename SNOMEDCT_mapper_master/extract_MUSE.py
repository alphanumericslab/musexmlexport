#!/usr/bin/env python

import numpy as np, os, os.path, sys
import pandas as pd
import itertools
from SNOMEDCT_mapper_master.mapper_muse import Mapper

#if __name__ == '__main__':

def extract_MUSE(dx_str):

    # Parse arguments.
#    if len(sys.argv) != 3:
#        raise Exception('Include the input and output directories as arguments, e.g., python driver.py input output.')

#    input_directory = sys.argv[1]
#    dfs = pd.read_excel(input_directory, sheet_name=None)
#    df = dfs['Sheet1']

    # iterating over rows using iterrows() function  
    nn=0

    file_id={}

    name_dx=pd.DataFrame([])

    data_total = []
    dx_total=[]
    dx_codes =list()


#    for i, j in df.iterrows():##itertools.islice(df.iterrows(), 10): 

#        name_dx.loc[i,'ID_file']=j['File_ID']
#    dx1=dx_str
    nn=0
    i=0

    dx_str=dx_str.split(",")
   
    for i,dx1 in enumerate(dx_str):
        dx1 = dx1.replace("[", "")
        dx1 = dx1.replace("]", "")

        if not pd.isnull(dx1) and not dx1==0:
            name_dx.loc[i,'Real_dx']=dx1
            snow_id=Mapper(dx1)
            df_snow=snow_id.standard_search()
            for k,m in df_snow.iterrows():
                nn+=1
                str_tmp="Term_"+str(nn)
                name_dx.loc[i,str_tmp]=m['Term']
                str_tmp="ID_"+str(nn)
                name_dx.loc[i,str_tmp]=m['ID']
                dx_codes.append(m['ID'])
    #            dx_code.append(m['Term'])
    #            dx_code.append(m['ID'])

    return dx_codes
    
#    name_dx.to_csv (r'Code_diagnosis_extracted.csv', index = False, header=True) #Don't forget to add '.csv' at the end of the path
