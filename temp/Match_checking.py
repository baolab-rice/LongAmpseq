import sys
import re

def main():
    with open(sys.argv[1],'r') as f:
        for line in f:
            cigar = line.split('\t')[5]
            x = re.search('M',cigar)
            y = re.findall("[0-9]+",cigar[:x.start()])
            print(y[-1])
    f.close()

if __name__ == '__main__':
    main()