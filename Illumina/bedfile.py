# python2.7
# ref: hg19

import sys
import math
import os
import numpy as np
import csv
import pandas as pd
import pickle
import argparse
import re
from sklearn import datasets
import matplotlib.pyplot as plt

def main():
    # use full path here
    bed_path = sys.argv[1]
    max_length = int(sys.argv[4])
    # bed_path = "/home/yp11/Documents/0219longamp/test_mergefiltered_2+.csv"

    df = pd.read_csv(bed_path, names=['chr', 'start', 'end', 'read_ID', 'score', 'strand'])

    # chr11	59136655	59136956	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
    # chr11	59138619	59138657	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
    # df.sort_values(by=['read_ID'])
    # ID_list = df.read_ID.unique()

    column_names = ['read_ID', 'start', 'length']
    df_all_filtered_output = pd.DataFrame(columns=column_names)
    df_all = pd.read_csv(bed_path, names=['chr', 'start', 'end', 'read_ID', 'score', 'strand'])
    id_all_list = df_all.read_ID.unique()

    # 3rd output: Include reads with more than 2 fragments, ignore reads with wrong directions
    for id in id_all_list:
        pos_list = []  # type: list
        chr_list = []
        dir_list = []
        pos_all_list = []
        chr_all_list = []
        dir_all_list = []
        df_id = df_all.loc[df_all['read_ID'] == id]
        for row in df_id.itertuples():
            chr_all_list.append(row.chr)
            pos_all_list.append(row.start)
            pos_all_list.append(row.end)
            dir_all_list.append(row.strand)
        # get the first and last
        if len(pos_all_list) >= 4:
            pos_list.append(pos_all_list[0])
            pos_list.append(pos_all_list[1])
            pos_list.append(pos_all_list[-2])
            pos_list.append(pos_all_list[-1])

            chr_list.append(chr_all_list[0])
            chr_list.append(chr_all_list[-1])

            dir_list.append(dir_all_list[0])
            dir_list.append(dir_all_list[-1])

        pos_list.sort()
        if chr_list[0] == chr_list[1] and 200 <= pos_list[2] - pos_list[1] < max_length and dir_list[0] == dir_list[1]:
            df_all_filtered_output = df_all_filtered_output.append(dict(
                read_ID=id, start=pos_list[1], length=pos_list[2] - pos_list[1]), ignore_index=True)
    df_all_filtered_output = df_all_filtered_output[df_all_filtered_output['length'] < max_length]
    df_all_filtered_output.to_csv(sys.argv[2], index=False)

    # consolidation_strict
    df_all_filtered_conso = df_all_filtered_output[['start', 'length']]
    df_all_filtered_conso_group = df_all_filtered_conso.groupby(df_all_filtered_conso.columns.tolist()).size().reset_index(). \
        rename(columns={0: 'repeat_num'})
    df_all_filtered_conso_group.to_csv(sys.argv[3], index=False)

    #Combine reads that are little off from the others
    #df_all_filtered_output

    for row in df_all_filtered_conso_group.itertuples():
        start = int(row.start)
        length = int(row.length)
        df_all_filtered_conso_group.at[row.Index, 'end'] = start + length


    # read numbers
    total_large_deletion_3 = len(df_all_filtered_output.index)
    # Open a file with access mode 'a'
    file_object = open('log.txt', 'a')
    file_object.write(bed_path)
    file_object.write('\n')
    file_object.write(str(total_large_deletion_3))
    file_object.write('\n')
    file_object.close()

if __name__ == '__main__':
    main()
