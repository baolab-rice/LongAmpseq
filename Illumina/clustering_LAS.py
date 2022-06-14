import sys

def read_file(filename):

    data_pair = []

    # Pay attention to the input format
    with open(filename,'r') as f:
        next(f)
        for line in f:
            line_split = line.split(',')
            start = int(line_split[1])
            deletion_size = int(line_split[2])
            data_pair.append([deletion_size, start])
    f.close()

    return data_pair

def cluster_samesize(listname,code):

    deletion_size_tolenrance = int(code.split('d')[0])
    deletion_start_tolenrance = int(code.split('d')[1].split('l')[0])

    list_sec = sorted(listname)

    list_filtered = []
    def find_cluster(listname,counter):

      list_filtered.append([[],[]])
      positions = [0]
      for i in range(len(listname)):
        if listname[i][0] - listname[0][0] <= deletion_size_tolenrance:
          if 0-deletion_start_tolenrance <= listname[i][1] - listname[0][1] <= deletion_start_tolenrance:
            list_filtered[counter][0].append(listname[i][0])
            list_filtered[counter][1].append(listname[i][1])
            positions.append(i)
      
      newlist = []
      for i in range(len(listname)):
        if i not in positions:
          newlist.append(listname[i])

      return newlist
    
    counter = 0
    while True:
      list_sec = find_cluster(list_sec,counter)
      counter += 1
      if len(list_sec) < 1:
        break

    list_final = []  
    for element in list_filtered:
      length = int(sum(element[0])/len(element[0]))
      start = int(sum(element[1])/len(element[1]))
      end = start + length
      cluster_size = len(element[1])
      list_final.append([cluster_size, start, end, length])

    return list_final

def output(listname):

    print("Cluster_size\tStart\tEnd\Length")
    for data in listname:
      print("{}\t{}\t{}".format(data[0],data[1],data[2]))

def main():

    dataset = read_file(sys.argv[1])
    dataset_list = cluster_samesize(dataset,"10d10l")
    output(dataset_list)

def write_output(listname, filename):

    f =  open(filename.replace(".txt", "_cluster.txt"),'w')
    f.write("Cluster_size\tStart\tEnd\tLength\n")  
    for data in listname:
      f.write("{}\t{}\t{}\t{}\n".format(data[0], data[1], data[2], data[3]))  
    f.close()

def cluster_generate(filename, code):

    dataset = read_file(filename)
    dataset_list = cluster_samesize(dataset, code)
    write_output(dataset_list, filename)

if __name__ == "__main__":
    main()

