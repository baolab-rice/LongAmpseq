# python2.7
# ref: hg19

import sys

def main():
    # read the output csv
    ld_path = sys.argv[1]
    
    f = open(ld_path)
    list_of_lines = f.readlines()[1:]
    f.close()

    ld_list = []
    for line in list_of_lines:
        ld_list.append([line.split(',')[0],line.split(',')[1],line.split(',')[2].strip()])

    # read the fastq
    fq_path = sys.argv[2]
    dict_fq = {}
    f = open(fq_path)
    list_of_lines = f.readlines()
    for i in range(len(list_of_lines)):
        if '>' in list_of_lines[i]:
            dict_fq[list_of_lines[i].split(' ')[0][1:]] = list_of_lines[i+1].strip()
    f.close()

    # mapping seq with read id
    for i in range(len(ld_list)):
        readid = ld_list[i][0]
        ld_list[i].append(dict_fq[readid])
    
    filename = open(ld_path,'w')
    filename.write('read_ID,start,length,seq')
    filename.write('\n')
    for read in ld_list:
        filename.write(read[0] + "," + read[1] + "," + read[2] + "," + read[3])
        filename.write('\n')

if __name__ == '__main__':
    main()
