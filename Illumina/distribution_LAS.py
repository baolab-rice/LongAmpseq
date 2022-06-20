#python3

#environment:
#conda install numpy scikit-learn pandas scipy
#conda install -c conda-forge matplotlib

import sys

import matplotlib.pyplot as plt
import pandas as pd

def red_blue_profile_pretty(df, cutsite, output_svg,linewidth,direct):
    df['start'] = df['start'].astype(float)
    df['length'] = df['length'].astype(float)

    df = df.sort_values(['start'], ascending=[0])
    df = df.sort_values(['length'], ascending=[1])

    df['start'] = df['start'] - cutsite

    i = 0

    # plt.rc('axes', titlesize=30) #change nothing
    plt.rc('axes', labelsize=20) # large label
    plt.rc('axes', linewidth=2)
    # plt.rc('lines', markersize=20) #change nothing
    plt.rc('xtick', labelsize=13)
    plt.tick_params(which='major', length=5, width=2, direction='inout')
    plt.yticks([])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    left = 0
    right = 0
    total = len(list(df.itertuples()))
    sim = 0

    for row in df.itertuples():
        mid = row.start + 0.5 * row.length
        end = row.start + row.length
        if direct == '-':
            temp1 = row.start
            temp2 = end
            nstart = temp2 * -1
            nend = temp1 * -1
            mid = mid * -1
        else:
            nstart = row.start
            nend = end
        if mid < 0:
            plt.hlines(i,nstart,nend, lw=linewidth, colors="red")
            left += 1
        else:
            plt.hlines(i,nstart,nend, lw=linewidth, colors="blue")
            right += 1
        i = i+1
        
    left_freq = str(left/total * 100)
    left_freq = left_freq.split('.')[0] + '.' + left_freq.split('.')[1][:2] + '%'
    right_freq = str(right/total * 100)
    right_freq = right_freq.split('.')[0] + '.' + right_freq.split('.')[1][:2] + '%'
    sim_rate = str(sim/total * 100)
    sim_rate = sim_rate.split('.')[0] + '.' + sim_rate.split('.')[1][:2] + '%'

    plt.xlabel("Relative position of large deletions")
    plt.xlim(-3000, 3000) # everthing except for HBG

    plt.text(-1500, 25, left_freq,
        verticalalignment='top', horizontalalignment='right',
        color='red', fontsize=20)

    plt.text(1500, 25, right_freq,
        verticalalignment='top', horizontalalignment='left',
        color='blue', fontsize=20)

    plt.text(1500, 5, "Sym:{}".format(sim_rate),
        verticalalignment='top', horizontalalignment='left',
        color='black', fontsize=15)

    plt.savefig(output_svg, bbox_inches='tight', format='svg')
    plt.clf()
    plt.close()

def read_LD(filename):
    ld_list = []
    with open(filename,'r') as f:
        next(f)
        for line in f:
            nline = line.split('\t')
            ld_list.append([nline[1], int(nline[2])-int(nline[1])])
        f.close()
    df = pd.DataFrame(ld_list, columns=['start', 'length'])

    return df

def distribution_generate(filename, cutsite):
    
    import matplotlib.pyplot as plt
    import pandas as pd

    df_group1 = read_LD(filename)
    output_svg = filename.replace(".txt", ".svg")
    red_blue_profile_pretty(df_group1, cutsite, output_svg, 2)

def main():

    df_group1 = read_LD(sys.argv[1])

    output_svg = sys.argv[1].replace(".csv", ".svg")

    cutsite = int(sys.argv[2])

    direct = sys.argv[3]

    red_blue_profile_pretty(df_group1, cutsite, output_svg, 2,direct)

if __name__ == '__main__':
    main()
