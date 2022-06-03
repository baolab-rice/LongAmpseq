import sys
import re

def main():
    with open(sys.argv[1],'r') as f:
        for line in f:
            cigar = line.split('\t')[5]
            x = re.split('[A-Z]',cigar)
            length = 0
            for num in x:
                if num != '':
                    length += int(num)
            start = int(line.split('\t')[3])
            if start + length > 2845:
                print(line.strip())
    f.close()

if __name__ == '__main__':
    main()