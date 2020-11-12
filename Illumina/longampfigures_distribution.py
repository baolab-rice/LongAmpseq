#python3
#ref: hg19

#environment:
#conda install numpy scikit-learn pandas scipy
#conda install -c conda-forge matplotlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import cycle
import math
import sys


# fuzzy consolidation with given width
# Assume input is sorted by starting position
# the row for combining, peak width, deletion length difference, the dataframe
# %%
# for greedy clustering
def combine_peaks(info_list, width, length, df_search):
    # info list format: dict
    # 5246101, 3845, 1
    # [['start', 'length', 'repeat_num']]

    width_temp_F = width
    length_temp_F = length
    width_temp_R = width
    length_temp_R = length

    while True:
        width_min = int(info_list.start) - width_temp_R
        width_max = int(info_list.start) + width_temp_F
        length_min = int(info_list.length) - length_temp_R
        length_max = int(info_list.length) + length_temp_F

        df_output = df_search[(df_search['start'].between(width_min, width_max)) &
                          (df_search['length'].between(length_min, length_max))]
        # get the maximum width and length
        width_max = df_output['start'].max()
        length_max = df_output['length'].max()
        width_min = df_output['start'].min()
        length_min = df_output['length'].min()

        df_expand = df_search[(df_search['start'].between(width_min - width, width_max + width)) &
                          (df_search['length'].between(length_min - length, length_max + length))]

        if (len(df_output.index) == len(df_expand.index)):
            break
        else:
            width_temp_F = width_max + width - int(info_list.start)
            length_temp_F = length_max + length - int(info_list.length)
            width_temp_R = int(info_list.start) - width_min + width
            length_temp_R = int(info_list.length) - length_min + length

    width_sum = 0
    length_sum = 0
    peak_total = 0
    for row in df_output.itertuples():
        width_sum = width_sum + int(row.start) * int(row.repeat_num)
        length_sum = length_sum + int(row.length) * int(row.repeat_num)
        peak_total = peak_total + int(row.repeat_num)
    peak_meanstart = width_sum / peak_total
    peak_meanlength = length_sum / peak_total

    peak_output = (peak_meanstart, peak_meanlength, peak_total)
    return df_output, peak_output


def delly_bed2group(df, chr, max_length):
    # cat ${file/.fastq/filtered_2+.bed} | tr "\\t" "," > ${file/.fastq/filtered_2+.csv}
    # assume the input is in csv
    # chr11	5247951	5248379	INV00000216
    column_names = ['start', 'length']
    df_group = pd.DataFrame(columns=column_names)
    df = df.loc[df['chr'] == chr]
    # df['type'] = df['type'].astype(str)
    df_del = df[df['type'].str.contains('DEL')]
    df_group['start'] = df_del['start']
    df_group['length'] = df_del['end'] - df_del['start']
    df_group = df_group.loc[df_group['length'] <= max_length]
    # definition of large deletion
    df_group = df_group.loc[df_group['length'] > 200]
    return df_group


def red_blue_profile(df, cutsite, output_eps):
    # df format: ['start', 'length']
    # df['start'] = df['start'].astype(float)
    # df['length'] = df['length'].astype(float)
    # df = df.sort_values(by=['start'])
    # print df
    # df = df.sort_values(by=['length'], ascending=False)
    # print df

    df = df.sort(['start'], ascending=[0])
    df = df.sort(['length'], ascending=[1])

    df['start'] = df['start'] - cutsite
    # print df

    i = 0

    for row in df.itertuples():
        mid = row.start + 0.5 * row.length
        # print 'mid:', mid
        end = row.start + row.length
        if mid < 0:
            plt.hlines(i, row.start, end, lw=0.1, colors="red")
        else:
            plt.hlines(i, row.start, end, lw=0.1, colors="blue")
        i = i+1

    plt.xlabel("Deletion Locations Spanning on Reference, 0 is cutsite")
    plt.ylabel("Large deletion patterns")
    plt.xlim(-4500, 4500)

    plt.savefig(output_eps, format='eps')


def red_blue_profile_pretty(df, cutsite, output_svg,linewidth):
    # df format: ['start', 'length']
    print('get to the function')
    df['start'] = df['start'].astype(float)
    df['length'] = df['length'].astype(float)
    # df = df.sort_values(by=['start'])
    # print df
    # df = df.sort_values(by=['length'], ascending=False)
    # print df

    df = df.sort_values(['start'], ascending=[0])
    df = df.sort_values(['length'], ascending=[1])

    df['start'] = df['start'] - cutsite
    # print df

    i = 0

    # plt.rc('axes', titlesize=30) #change nothing
    plt.rc('axes', labelsize=20) # large label
    plt.rc('axes', linewidth=2)
    # plt.rc('lines', markersize=20) #change nothing
    plt.rc('xtick', labelsize=15)
    plt.tick_params(which='major', length=5, width=2, direction='inout')
    plt.yticks([])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    for row in df.itertuples():
        # print(row)
        mid = row.start + 0.5 * row.length
        # print 'mid:', mid
        end = row.start + row.length
        if mid < 0:
            plt.hlines(i, row.start, end, lw=linewidth, colors="red")
        else:
            plt.hlines(i, row.start, end, lw=linewidth, colors="blue")
        i = i+1

    plt.xlabel("Relative position of large deletions")
    plt.xlim(-3000, 3000) # everthing except for HBG
    # plt.xlim(-8000, 2000)

    plt.savefig(output_svg, bbox_inches='tight', format='svg')

def red_blue_profile_cellline(df, cutsite, output_eps):
    # df format: ['start', 'length']
    # df['start'] = df['start'].astype(float)
    # df['length'] = df['length'].astype(float)
    # df = df.sort_values(by=['start'])
    # print df
    # df = df.sort_values(by=['length'], ascending=False)
    # print df

    df = df.sort(['start'], ascending=[0])
    df = df.sort(['length'], ascending=[1])

    # for cell line we need to flip the orientation
    # df['start'] = - df['start'] + cutsite
    df['start'] = df['start'] - cutsite
    # print df

    i = 0

    for row in df.itertuples():
        mid = row.start + 0.5 * row.length
        # print 'mid:', mid
        end = row.start + row.length
        if mid > 0:
            plt.hlines(i, -end, -row.start, lw=0.1, colors="red")
        else:
            plt.hlines(i, -end, -row.start, lw=0.1, colors="blue")
        i = i+1

    plt.xlabel("Deletion Locations Spanning on Reference, 0 is cutsite")
    plt.ylabel("Large deletion patterns")
    plt.xlim(-8500, 2500)

    plt.savefig(output_eps, format='eps')


def v_distribution(samplename, type, df1, df2, cutsite):
    df1['start'] = df1['start'] - cutsite
    df1['mid'] = df1['start'] + 0.5 * df1['length']
    label1 = samplename + ' ' + type + '1'

    df2['start'] = df2['start'] - cutsite
    df2['mid'] = df2['start'] + 0.5 * df2['length']
    label2 = samplename + ' ' + type + '2'

    alpha = 0.5

    plt.scatter(x=df1['mid'], y=df1['length'], s=1, label=label1, marker='o', color=[1., alpha, alpha])
    plt.scatter(x=df2['mid'], y=df2['length'], s=0.05, label=label2, marker='x', color=[alpha, alpha, 1.])
    plt.xlabel("Deletion Middle Point Locations")
    plt.ylabel("Deletion Lengths")
    plt.legend(fontsize = "xx-small")
    plt.xlim(-2000, 2000)
    plt.ylim(0, 4000)

    output_eps = samplename + ' ' + type + '.eps'
    plt.savefig(output_eps, format='eps')
    plt.clf()


def v_distribution_glob(samplename, label1, label2, df1, df2, cutsite):
    df1['start'] = df1['start'] - cutsite
    df1['mid'] = df1['start'] + 0.5 * df1['length']

    df2['start'] = df2['start'] - cutsite
    df2['mid'] = df2['start'] + 0.5 * df2['length']

    alpha = 0.5

    plt.scatter(x=df1['mid'], y=df1['length'], s=1, label=label1, marker='o', color=[1., alpha, alpha])
    plt.scatter(x=df2['mid'], y=df2['length'], s=0.05, label=label2, marker='x', color=[alpha, alpha, 1.])
    plt.xlabel("Deletion Middle Point Locations")
    plt.ylabel("Deletion Lengths")
    plt.legend(fontsize = "xx-small")
    plt.xlim(-2000, 2000)
    plt.ylim(0, 4000)

    output_eps = samplename + '.eps'
    plt.savefig(output_eps, format='eps')
    plt.clf()


def v_distribution_glob_nanoill(samplename, label1, label2, df1, df2, cutsite):
    df1['start'] = df1['start'] - cutsite
    df1['mid'] = df1['start'] + 0.5 * df1['length']

    df2['start'] = df2['start'] - cutsite
    df2['mid'] = df2['start'] + 0.5 * df2['length']

    alpha = 0.5

    plt.scatter(x=df1['mid'], y=df1['length'], s=0.1, label=label1,alpha=0.3)
    plt.scatter(x=df2['mid'], y=df2['length'], s=0.1, label=label2,alpha=0.3)
    plt.xlabel("Deletion Middle Point Locations")
    plt.ylabel("Deletion Lengths")
    plt.legend(fontsize = "xx-small")
    plt.xlim(-2000, 2000)
    plt.ylim(0, 5300)

    output_eps = samplename + '.eps'
    plt.savefig(output_eps, format='eps')
    plt.clf()


def del_to_group(df_del,output):
    # format of df_del
    # df1 = pd.read_csv('NoiseCanceledNanopore.csv', names=['start', 'length'])
    # consolidation_strict
    df_group = df_del.groupby(df_del.columns.tolist()).size().reset_index(). \
        rename(columns={0: 'repeat_num'})
    df_group.to_csv(output, index=False)



def main():

    input = sys.argv[1]
    output_eps = input.replace(".csv", ".eps")
    output_svg = input.replace(".csv", ".svg")

    if 'BCL11A' in input:
        cutsite = 60722401
        length = 4267
        chr = 'chr2'
    elif 'C19' in input:
        cutsite = 1004
        length = 9423
        chr = 'GFPBFP'
    elif 'OT18' in input:
        cutsite = 59137387
        length = 5045
        chr = 'chr11'
    elif 'R66' in input or 'HBB' in input:
        cutsite = 5248229
        length = 5465
        chr = 'chr11'
    elif 'R02' in input:
        cutsite = 5248214
        length = 5465
        chr = 'chr11'
    elif 'HBG1' in input:
        cutsite = 2759
        length = 6718
        chr = 'HBGLONGAMP'
    elif len(sys.argv) ==4:
        chr = str(sys.argv[2])
        cutsite = int(sys.argv[3])
        length = int(sys.argv[4])
    else:
        cutsite = 0
        length = 0
        chr = ''

    # just to get read numbers before clustering
    df1 = pd.read_csv(input)
    # df1 = df1.loc[df1['repeat_num'] > 1]
    print(input)
    print(len(df1.index))

# HBG need to be +-4500

    # # /////////////////////////////////////////////////////////
    # df1 = pd.read_csv(input)
    # df1 = df1.loc[df1['repeat_num'] > 1]
    # red_blue_profile_pretty(df1, cutsite, output_svg, 1)
    # # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    #optimize figure drawing
    # input = 'R66SWT-2_S6_L001_R1_001_merged30_filtered_sorted.csv'
    df_bed1 = pd.read_csv(input, names=['chr', 'start', 'end', 'type'])
    # df_bed1 =
    # pd.read_csv('R66SWT-2_S6_L001_R1_001_merged30_filtered_sorted.csv',  names=['chr', 'start', 'end', 'type'])
    df_group1 = delly_bed2group(df_bed1, chr, length)
    red_blue_profile_pretty(df_group1, cutsite, output_svg, 2)

    # /////////////////////////////////////////////////////////


    # /////////////////////////////////////////////////////////
    # compare R02 and R66 result
    # R66SHIFI-2_largelargedel_group.csv
    # R66SHIFI-1_501largedel_group.csv ( likely mislabeled)
    # df1 = pd.read_csv('R66SHIFI-2_501largedel_group.csv')
    # df1 = df1.loc[df1['repeat_num'] > 1]
    # df2 = pd.read_csv('R02-2_501largedel_group.csv')
    # df2 = df2.loc[df2['repeat_num'] > 1]
    # df3 = pd.read_csv('R66SHIFI-1_501largedel_group.csv')
    # df3 = df3.loc[df3['repeat_num'] > 1]
    #v_distribution_glob('R66S2vsR02', 'R02', 'R66S2', df2, df1, cutsite)
    # v_distribution_glob('R66S2vsR66S1', 'R66S1', 'R66S2', df3, df1, cutsite)

    # /////////////////////////////////////////////////////////
    # redraw dbscan result after noise cancellation
    # data from Yilei
    # df1 = pd.read_csv('NoiseCanceledNanopore.csv', names=['start', 'length'])
    # df1['start'] = df1['start'] + cutsite
    # df1_group = df1.groupby(df1.columns.tolist()).size().reset_index(). \
    #     rename(columns={0: 'repeat_num'})
    #
    # df2 = pd.read_csv('NoiseCanceledIllumina.csv', names=['start', 'length'])
    # df2['start'] = df2['start'] + cutsite
    # df2_group = df2.groupby(df2.columns.tolist()).size().reset_index(). \
    #     rename(columns={0: 'repeat_num'})
    #
    # df1_group.to_csv(r'NoiseCanceledNanopore_dedup.csv')
    # df2_group.to_csv(r'NoiseCanceledIllumina_dedup.csv')
    #
    # v_distribution_glob_nanoill('R66SHIFI-2_DBSCAN', 'Nanopore Sequenced Sample', 'Illumina Sequenced Sample', df1_group, df2_group, cutsite)
    # s1 = pd.merge(df1_group, df2_group, how='inner', on=['start', 'length'])
    # s1.to_csv(r'Nanoill_overlap.csv')

    # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    # compare small and large run reads visualization:
    # R66SHIFI-2_largelargedel_group.csv
    # R66SHIFI-1_501largedel_group.csv ( likely mislabeled)
    # df1 = pd.read_csv('R66SHIFI-2_largelargedel_group.csv')
    # df1 = df1.loc[df1['repeat_num'] > 1]
    # df2 = pd.read_csv('R66SHIFI-1_501largedel_group.csv')
    # df2 = df2.loc[df2['repeat_num'] > 1]
    # v_distribution_glob('R66SHIFI-2', 'Low coverage', 'High coverage', df2, df1, cutsite)
    # s1 = pd.merge(df1, df2, how='inner', on=['start', 'length'])
    # s1.to_csv(r'nocluster_overlap.csv')

    # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    # show delly comparison of small and large illumina run
    # pwd = ~/Documents/0219longamp/0309_fig4C/0315

    # df_bed1 = pd.read_csv(input,  names=['chr', 'start', 'end', 'type'])
    # df_group1 = delly_bed2group(df_bed1, chr, length)
    # input2 = 'R66SHIFI-1_S3_L001_R1_001_merged30_filtered_sorted.csv'
    # df_bed2 = pd.read_csv(input2, names=['chr', 'start', 'end', 'type'])
    # df_group2 = delly_bed2group(df_bed2, chr, length)
    # v_distribution_glob('R66SHIFI_delly_0427', 'Low coverage', 'High coverage', df_group2, df_group1, cutsite)
    # #
    # s1 = pd.merge(df_bed1, df_bed2, how='inner', on=['chr', 'start', 'end'])
    # s1.to_csv(r'delly_overlap.csv')
    # /////////////////////////////////////////////////////////


    # /////////////////////////////////////////////////////////
    # V comparisons
    # sample1_502 = input.replace("501", "502")
    #
    # sample2_501 = input.replace("-1_", "-2_")
    # print sample2_501
    # sample2_502 = sample1_502.replace("-1_", "-2_")
    # print sample2_502
    # df_1_501 = pd.read_csv(input)
    # df_1_501 = df_1_501.loc[df_1_501['repeat_num'] >1]
    #
    # df_1_502 = pd.read_csv(sample1_502)
    # df_1_502 = df_1_502.loc[df_1_502['repeat_num'] >1]
    #
    # df_2_501 = pd.read_csv(sample2_501)
    # df_2_501 = df_2_501.loc[df_2_501['repeat_num'] >1]
    #
    # df_2_502 = pd.read_csv(sample2_502)
    # df_2_502 = df_2_502.loc[df_2_502['repeat_num'] >1]
    #
    # # strangely tech rep and bio rep cannot run together
    # # type = 'Technical replicate'
    # # samplename1 = input.replace("-1_501largedel_group.csv", "_1")
    # # samplename2 = input.replace("-1_501largedel_group.csv", "_2")
    # # v_distribution(samplename1, type, df_1_501, df_1_502, cutsite)
    # # v_distribution(samplename2, type, df_2_501, df_2_502, cutsite)
    #
    # type = 'Biological replicate'
    # samplename1 = input.replace("-1_501largedel_group.csv", "_501")
    # samplename2 = input.replace("-1_501largedel_group.csv", "_502")
    # v_distribution(samplename1, type, df_1_501, df_2_501, cutsite)
    # v_distribution(samplename2, type, df_1_502, df_2_502, cutsite)
    #
    #

    # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    # draw cell line data
    # df_group = pd.read_csv(input)
    # df_group = df_group.loc[df_group['repeat_num'] >1]
    # red_blue_profile_cellline(df_group, cutsite, output_eps)

    # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    # just draw standard profiles
    # df_group = pd.read_csv(input)
    # df_group = df_group.loc[df_group['repeat_num'] >1]
    # red_blue_profile(df_group, cutsite, output_eps)

    # /////////////////////////////////////////////////////////

    # /////////////////////////////////////////////////////////
    # for file in *_sorted.csv; do python ~/Documents/scripts/0314_longampfigures_distribution.py ${file}; done
    # draw delly profiles.


    # df_bed = pd.read_csv(input,  names=['chr', 'start', 'end', 'type'])
    #
    # df_group = delly_bed2group(df_bed, chr, length)
    # red_blue_profile(df_group, cutsite, output_eps)
    # /////////////////////////////////////////////////////////



if __name__ == '__main__':
    main()

