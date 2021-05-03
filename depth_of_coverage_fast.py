# Reads depth by bp from a text file
# Outputs a .csv containing sequence name, breadth of coverage, sum of depths,
# length, depth of coverage, and bases mapped.
# Written by Zachary A. Batz

import sys
import getopt
import csv
from timeit import default_timer as timer
import re

def usage():
    print('Usage: python depth_of_coverage_fast.py -t <textfile.txt> -o <output.csv>')

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "ht:o:", \
            ["help", "text=", "output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    txtfile = None
    outputfile = None
    table = False

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-t", "--txt"):
            txtfile = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
        else:
            assert False, "unhandled option"

    if (txtfile == None) or (outputfile == None):
        usage()
        sys.exit()

    try:
        start = timer()
        ## Key = gene_exon
        ## Value = [length, positions_mapped, total_depth]
        reads_dict = {}
        genome_data = [0,0,0]
        with open(txtfile) as inf:
            reader = csv.reader(inf, delimiter = "\t")
            for bp in reader:
                key = bp[0]

                if key not in reads_dict.keys():
                    reads_dict[key]=[0,0,0]

                value = reads_dict[key]


                ## Add position to length
                value[0] += 1
                genome_data[0] += 1

                ## If mapped, add to positions mapped
                if int(bp[2]) > 0:
                    value[1] += 1
                    genome_data[1] += 1

                ## Add depth to total depth
                value[2] += int(bp[2])
                genome_data[2] += int(bp[2])

        ## sort keys alphabetically
        sorted_reads_dict = dict( sorted(reads_dict.items(), key=lambda x: x[0].lower()) )

        with open(outputfile,'w') as ofile:
            writer = csv.writer(ofile)

            header = ["Sequence Name","Length","Bases Mapped","Breadth of Coverage","Sum of Depths","Depth of Coverage"]
            writer.writerow(header)

            for key, value in sorted_reads_dict.items():
                length = value[0]
                bmapped = value[1]
                broc = round((float(bmapped)/length)*100,10)
                sod = value[2]
                doc = round((float(sod)/length)*100,9)

                out_line = [key,length,bmapped,broc,sod,doc]
                for i in range(0,len(out_line)):
                    if out_line[i] == 0:
                        out_line[i] = "NA"
                writer.writerow(out_line)


    except Exception as e:
        print("An exception occurred:")
        raise
        sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])
