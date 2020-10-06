# Reads sequences from a .bam
# Outputs a .csv containing sequence name, breadth of coverage, sum of depths,
# length, depth of coverage, and bases mapped. Can use a '-t' flag to output in
# a friendly Excel table format
#
# Author: Matt Preston
# Created on: Jan 18, 2016

import sys
import getopt
import depth_of_coverage # coverage

def usage():
    print 'Usage: python depth_of_coverage_impl.py -b <bamfile.bam> ' \
                 '-o <output.csv> [-t]'

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hb:o:t", \
            ["help", "bam=", "output=", "table"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    bamfile = None
    outputfile = None
    table = False
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-b", "--bam"):
            bamfile = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ("-t", "--table"):
            table = True
        else:
            assert False, "unhandled option"
    
    if (bamfile == None) or (outputfile == None):
        usage()
        sys.exit()
    
    try:
        # Get list of tuples from .bam
        # i.e. (maybe not in this order however)
        # [(SeqName, "Breadth of Coverage", b.o.c.),
        #  (SeqName, "Sum of Depths", s.o.d.),
        #  (SeqName, "Length", len),
        #  (SeqName, "Depth of Coverage", d.o.c.),
        #  (SeqName, "Bases Mapped", b.m.),
        #  (SeqName2, "Breadth of Coverage", b.o.c.2), ...]
        coverage = depth_of_coverage.coverage(bamfile)
        order = ["Length","Bases Mapped","Breadth of Coverage",
                 "Sum of Depths","Depth of Coverage"]
        # Sort coverage by the order of colums desired
        coverage = [entry for q in range(0,len(coverage),5)  \
                            for i in order                   \
                                for entry in coverage[q:q+5] \
                                    if entry[1] == i]
        with open(outputfile, 'w') as f:
            if table: # Excel table format
                f.write("Sequence Name,")
                for i in order:
                    f.write("%s," % i)
                f.write("\r\n")
                # Compress more [(SeqName,len,b.m.,b.o.c.,s.o.d.,d.o.c.)...]
                templist = []
                for q in range(0, len(coverage),5):
                    if coverage[q][0] != "genome":
                        tup = (coverage[q][0],)
                        # If i[2] == 0: return "NA" else return i[2]
                        tup += tuple(i[2] if i[2] else "NA" \
                                     for i in coverage[q:q+5])
                        templist.append(tup)
                # Sort by SeqName before writing
                templist.sort(key=lambda entry: entry[0].lower())
                for entry in templist:
                    f.write("%s,%s,%s,%s,%s,%s\r\n" % entry)
            else: # Literal text
                for i in coverage:
                    f.write("%s,%s,%s\r\n" % i)
    except Exception as e:
        print "An exception occurred:"
        raise
        sys.exit()
        
if __name__ == '__main__':
    main(sys.argv[1:])
        
