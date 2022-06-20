import sys

def read_from_LIhits(filename):
    hits = {}
    distance = 3000
    Start = 0
    Chr_marker = 'Chr'
    with open(filename,'r') as f:
        next(f)
        for line in f:
            splitline = line.strip().split('\t')
            Chr = splitline[1]
            if Chr == Chr_marker:
                nStart = int(splitline[2])
                if Start == 0:
                    Start = nStart
                    ID = splitline[0] 
                    hits[ID] = (Chr,Start)
                # For now just keep the first hit
                elif nStart - Start > distance:
                    Start = 0
            else:
                Chr_marker = Chr
                Start = 0
    f.close()

    return hits

def output(dictname,filename):
    
    outputfile = filename[:-4] + '_binned.txt'
    f = open(outputfile,'w')
    for ID in dictname.keys():
        f.write(ID)
        f.write('\t')
        f.write(dictname[ID][0])
        f.write('\t')
        f.write(str(dictname[ID][1]))
        f.write('\n')
    f.close()

def main():
    filename = sys.argv[1]
    hits = read_from_LIhits(filename)
    output(hits,filename)

if __name__ == '__main__':
    main()
    
