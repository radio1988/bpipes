from fuzzywuzzy import fuzz
from itertools import zip_longest as zip
import gzip
import sys
import argparse
import re


def parse_args(argv):
    parser = argparse.ArgumentParser(description = 'Help information')
    parser.add_argument('file1', help='xxx.R1.fq')
    parser.add_argument('file2', help='xxx.R2.fq') 
    parser.add_argument('-s', '--min_similarity', help="matched_length/full_length", 
                        default=61.5, type=float)
    
    return parser.parse_args(argv)


# Note:
print('''
>> Desc:
- Keep reads if they have 26bp adaptors in the beginning of R2. Allowing similarity as low as min_similarity (default 61.5)
- Trim the first 26bp in R2 of the qualified reads
- In PE mode 
''')

# Parameters
args = parse_args(sys.argv[1:])
file1 = args.file1
file2 = args.file2

out1=file1.replace(".fastq.gz", "")
out1=out1+".flt.fastq.gz"

out2=file2.replace(".fastq.gz", "")
out2=out2+".flt.fastq.gz"

adaptor="GTTCAGAGTTCTACAGTCCGACGATC"
min_similarity=float(args.min_similarity) #61.5  # 16/26

print("inputs:", file1, file2)
print("similarity", min_similarity)
print("outputs:", out1, out2)

# Prep
out_x = [1,2,3,4]
out_y = [1,2,3,4]

write_x = gzip.open(out1, 'wt')
write_y = gzip.open(out2, "wt")

# Work
# Work
with gzip.open(file1, "rt") as f1, gzip.open(file2, "rt") as f2: #rt for text mode
    i=0
    o=0
    for x, y in zip(f1, f2):  # read two files parrallel 
        i += 1
        
        # residue 0,1,2,3
        residue = i % 4
        if residue == 0:
            residue = 4
        residue = residue -1      
        
        # clean up
        x = x.strip()
        y = y.strip()
        y_head = y[:26]
        
        # store entry for one read
        out_x[residue] = x
        out_y[residue] = y
    
        
        # check R1, R2 match
        if  residue == 0:  # ID
            if not x.split()[0] == y.split()[0]:
                raise Exception ("R1 R2 fastq id not matching, maybe misaligned R1 R2 files")
                
        # check read content
        if residue == 1:  # SEQ
            # check adaptor in R2
            similarity = fuzz.ratio(y_head, adaptor)
            adap_qual = ("Good" if similarity > min_similarity else "Bad")
            
            # check tailing As in R2: 2-39
            p_taila = re.compile('A+$')
            if p_taila.search(y):
                len_taila = len(p_taila.search(y).group())
                tailing_as = ("Good" if (len_taila >1 and len_taila < 40) else "Bad")
            else:
                len_taila = 0
                tailing_as = "Bad"
            
            # check leading with more than one T in R1
            p_leading_ts = re.compile('^T+')
            if p_leading_ts.search(x):
                len_leading_ts=len(p_leading_ts.search(x).group())
                leading_ts="Good"
                # print(x, len_leading_ts)
            else:
                len_leading_ts=0
                leading_ts="Bad"


        # Output
        if residue == 3 and adap_qual=="Good" and leading_ts=="Good" and tailing_as=="Good":
#             print ("read number", i//4, "\t similarity:{}\tlen_R1_ts:{}\tlen_R2_as:{}\n".
#                    format(similarity,len_leading_ts, len_taila))
            
            # trim R2 adaptor
            out_y[1] = out_y[1][26:]
            out_y[3] = out_y[3][26:]
            # trim R2 tailing As
            out_y[1] = out_y[1][:len(out_y[1])-len_taila]
            out_y[3] = out_y[3][:len(out_y[3])-len_taila]
            # trim R1 leading Ts
            out_x[1] = out_x[1][len_leading_ts:]
            out_x[3] = out_x[3][len_leading_ts:]
            
            # next if nothing left
            if (len(out_x[1])<13 or len(out_y[1])<13):
#                 print ("too short after trimming")
                next

            
#             # output preview
#             print("R1:")
#             print("\n".join(out_x))
#             print("\n")
#             print("R2")
#             print("\n".join(out_y))
#             print("\n")
            
            # output
            write_x.write("\n".join(out_x))
            write_x.write("\n")
            write_y.write("\n".join(out_y))
            write_y.write("\n")
            o += 1
            
        
#         if i > 10: 
#             raise Exception ("too long")
#             #print (i//4, fuzz.ratio(y.strip(), adaptor), y.strip())
#         #print("{0}\t{1}\t{2}".format(i, x.strip(), y.strip()))




write_x.close()
write_y.close()

print(">>>Finished, there are {} read_pairs passed filter, trimming".format(o))