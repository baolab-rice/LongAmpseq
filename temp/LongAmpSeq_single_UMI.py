from subprocess import Popen, PIPE
import sys 
import argparse
import re

### /Users/mingmingcao/Desktop/LongAmpSeq/Single_UMI/5_9_22_synthetic_library\ files/aligned\ bam\ sam\ files/pcr3_6x_r1_sorted_filtered_length_region.sam

"""
def get_umi_bin_list(_dir):

    arguments = ['ls {}/raconx2/ | xargs'.format(_dir)]
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)

    umibins_raw = stdout.decode("utf-8")
    umibins_names = [x.strip() for x in umibins_raw.split(" ") if "umi" in x]
    umibins_dirs = [_dir + "/raconx2/" + x.strip() for x in umibins_raw.split(" ") if "umi" in x]

    return umibins_names,umibins_dirs
"""
## Raw read alignment

## Read trimming and filtering 


### UMI & Barcode extraction 
### UMI filtering 
### Barcode filtering 
### UMI & Barcode combanition 
### Raw read alignment 
### UMI binning 
### Consensus seq generation
### Coverage filtering (singletons)
global BC_list

BC_list = ["ACCATCACG",
           "ACCCGATGT",
           "ACCTTAGGC",
           "ACCTGACCA",
           "ACCACAGTG",
           "ACCGCCAAT",
           "ACCCAGATC",
           "ACCACTTGA",
           "ACCGATCAG",]

def readfile(filename):
    reads_dict = {}
    counter = 0
    BC_L = "CAGGAAAC"
    counter = 0
    with open(filename,'r') as f:
        for line in f:
            if line.startswith("@"):
                pass
            else:
                title = line.split('\t')[0]
                # For different references, the start postions are different 
                # start = 41 - int(line.split('\t')[3])
                # read = line.split('\t')[9][start:]
                x = re.search(BC_L, line.split('\t')[9])
                if x != None and x.start() > 17:
                    read = line.split('\t')[9][x.start()-18:]
                # Filter for length from UMI to BC
                    if len(read) > 118:
                        reads_dict[title] = read

                y= re.search("CACCTGCT",line.split('\t')[9])
                if y != None:
                    counter += 1
                
                
    print(counter)

    f.close()

    return reads_dict

def UMI_and_BC_searching(read):   
    regex = "[ATCG]{3}[CT]{1}[AG]{1}[ATCG]{3}[CT]{1}[AG]{1}[ATCG]{3}[CT]{1}[AG]{1}[ATCG]{3}"
    x = re.search(regex,read)
    return x


def main():
    reads = readfile(sys.argv[1])
    del_read_list = []
    frag1 = {}
    frag2 = {}
    frag3 = {}
    frag4 = {}
    frag5 = {}
    frag6 = {}
    frag7 = {}
    frag8 = {}
    frag9 = {}
    for title, read in reads.items():
        x = UMI_and_BC_searching(read[:18]) 
        if x == None or x.start() != 0:
            del_read_list.append(title)
        else:
            for i in range(len(BC_list)):
                y = re.search(BC_list[i],read)
                if y != None:
                    if y.start() - x.start() > 50: # Some BCs are near UMI
                        if i == 0:
                            frag1[title] = read[x.start():x.end()]
                        if i == 1:
                            frag2[title] = read[x.start():x.end()]
                        if i == 2:
                            frag3[title] = read[x.start():x.end()]
                        if i == 3:
                            frag4[title] = read[x.start():x.end()]
                        if i == 4:
                            frag5[title] = read[x.start():x.end()]
                        if i == 5:
                            frag6[title] = read[x.start():x.end()]
                        if i == 6:
                            frag7[title] = read[x.start():x.end()]
                        if i == 8:
                            frag9[title] = read[x.start():x.end()]
                else:
                    if "CACCTGCT" in read:
                        frag8[title] = read[x.start():x.end()]
    for read in del_read_list:
        reads.pop(read, None)
    
    # Write to the file
    Frag_List = [frag1,frag2,frag3,frag4,frag5,frag6,frag7,frag8,frag9,]
    count = 0
    for i in range(9):
        with open("Fragment{}_UMIs.txt".format(i+1), "w") as f:
            for title,read in Frag_List[i].items():
                f.write(">" + title)
                f.write('\n')
                f.write(read)
                f.write('\n')
        print(len(Frag_List[i]))
        count = count + len(Frag_List[i])
    print("The total UMIs is:")
    print(count)


if __name__ == '__main__':
    main()


#1. Applied strict pattern for NNNYR
