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
    parser.add_argument('-t', '--trim', help="to trim leading ts in R1 and tailing as in R2 or not {0,1}, default 1", 
                    default=1, type=int)
    
    return parser.parse_args(argv)


def search_leading_ts(seq):
    p = re.compile('^T+')
    if p.search(seq):
        length=len(p.search(x).group())
        status="Good"
        # print(x, len_leading_ts)
    else:
        length=0
        status="Bad"

    return(status, length)


def search_tailing_as(seq):
    p = re.compile('A+$')
    if p.search(seq):
        length = len(p.search(seq).group())
        status = ("Good" if (length >1 and length < 40) else "Bad")
    else:
        length = 0
        status = "Bad"

    return(status, length)



# Note:
print('''
<find_trim_3P_adaptor.py> PE mode
Help Information:

Example: 
python3 find_trim_3P_adaptor.py R1 R2 -s 61.5 -t 1

>> Workflow
> Read R1 and R2 in PE mode

> Filtering
1. Contains adaptor in the beginning of R2: If adaptor similarity > 61.5% for the beginning 26 bp in R2
2. R2 ends with 2-39 terminal A's 
2. R1 starts with leading T's in the beginning

> Trimming
1. removed first 26bp in R2 (always)
2. removed ending As in R2 (optional by -t 1)
3. removed starting Ts in R1 (optinoal by -t 1)

>> Update Note (new to old):
- added progress info
- remove pairs that are too short after trimming
- print args

>> Todo: 
1. More flexible

''')




# Parameters
args = parse_args(sys.argv[1:])

out1=args.file1.replace(".fastq.gz", "")
out1=out1+".flt.fastq.gz"

out2=args.file2.replace(".fastq.gz", "")
out2=out2+".flt.fastq.gz"

adaptor="GTTCAGAGTTCTACAGTCCGACGATC"
min_similarity=float(args.min_similarity) #61.5  # 16/26

print(">>> Began:")
print("cmd: ", args)
print("inputs:", args.file1, args.file2)
print("similarity", min_similarity)
print("outputs:", out1, out2)

# Prep
out_x = [1,2,3,4]
out_y = [1,2,3,4]

write_x = gzip.open(out1, 'wt')
write_y = gzip.open(out2, "wt")

# Work
with gzip.open(args.file1, "rt") as f1, gzip.open(args.file2, "rt") as f2: #rt for text mode
    i=0
    o=0
    for x, y in zip(f1, f2):  # read two files parrallel 
        i += 1
        
        # residue 0,1,2,3
        residue = i % 4
        if residue == 0:
            residue = 4
        residue = residue - 1      

        # progress report
        if residue == 4 and i//4 % 1000000 == 1: 
            print("{} M pairs processed".format(i//4000000))
        
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
            (tailing_as, len_taila) = search_tailing_as(y)
            
            # check leading with more than one T in R1
            (leading_ts, len_leading_ts) = search_leading_ts(x)


        # Output at line 4
        if residue == 3 and adap_qual=="Good" and leading_ts=="Good" and tailing_as=="Good":
            # print ("read number", i//4, "\t similarity:{}\tlen_R1_ts:{}\tlen_R2_as:{}\n".
            #        format(similarity,len_leading_ts, len_taila))
            
            # trim R2 adaptor
            out_y[1] = out_y[1][26:]
            out_y[3] = out_y[3][26:]

            if args.trim:
                # trim R2 tailing As
                out_y[1] = out_y[1][:len(out_y[1])-len_taila]
                out_y[3] = out_y[3][:len(out_y[3])-len_taila]
                # trim R1 leading Ts
                out_x[1] = out_x[1][len_leading_ts:]
                out_x[3] = out_x[3][len_leading_ts:]
            
            # if not too short after trimming
            if not (len(out_x[1])<13 or len(out_y[1])<13):
                # output
                write_x.write("\n".join(out_x))
                write_x.write("\n")
                write_y.write("\n".join(out_y))
                write_y.write("\n")
                o += 1
                # output preview
                # print("R1:")
                # print("\n".join(out_x))
                # print("\n")
                # print("R2")
                # print("\n".join(out_y))
                # print("\n")




print(">>>Finished, there are {}/{} read_pairs passed".format(o, i//4))
