#python3

#environment:
#conda install numpy scikit-learn pandas scipy
#conda install -c conda-forge matplotlib

import sys

import matplotlib.pyplot as plt
import pandas as pd

def red_blue_profile_pretty(df, cutsite, output_svg,linewidth,direction):
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

    for row in df.itertuples():
        mid = row.start + 0.5 * row.length
        end = row.start + row.length
            if direction == '-':
                row.start = temp1
                end = temp2
                row.start = -temp2
                end = -temp1
        if mid < 0:
            plt.hlines(i, row.start, end, lw=linewidth, colors="red")
        else:
            plt.hlines(i, row.start, end, lw=linewidth, colors="blue")
        i = i+1

    plt.xlabel("Relative position of large deletions")
    plt.xlim(-3000, 3000) # everthing except for HBG

    plt.savefig(output_svg, bbox_inches='tight', format='svg')
    plt.clf()
    plt.close()

def read_LD(filename):
    ld_list = []
    with open(filename,'r') as f:
        next(f)
        for line in f:
            nline = line.split(',')
            ld_list.append([nline[1], int(nline[1]) + int(nline[2]), nline[2]])
        f.close()
    df = pd.DataFrame(ld_list, columns=['start', 'end', 'length'])

    return df

def distribution_generate(filename, cutsite):
    
    import matplotlib.pyplot as plt
    import pandas as pd

    df_group1 = read_LD(filename)
    output_svg = filename.replace(".txt", ".svg")
    red_blue_profile_pretty(df_group1, cutsite, output_svg, 2)

def main():

    df_group1 = read_LD(sys.argv[1])

    output_svg = sys.argv[1].replace(".txt", ".svg")

    cutsite = int(sys.argv[2])

    direction = sys.argv[3]

    red_blue_profile_pretty(df_group1, cutsite, output_svg, 2,direction)

    

if __name__ == '__main__':
    main()
