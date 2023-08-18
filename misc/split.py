import gzip,sys,os

CODES = {
'100KME':'GGACTCCT+AGCAATCC',
'100KAC':'GGACTCCT+TCCTTGAC',
'500KME':'TAGGCATG+AGCAATCC',
'500KAC':'TAGGCATG+TCCTTGAC',
'1xK4':'CGTACTAG+AGGCTTAC',
'2xK4':'AGGCAGAA+AGGCTTAG',
'4xK4':'TCCTGAGC+AGGCTTAG',
}

test = False


def barcode_dist (barcode, code):
    '''
    NAGNCTTN+TCCTTGAC
    '''
    score = 0
    length = 0
    for i,s in enumerate(barcode):
        b = barcode[i]
        c = code[i]
        if b == 'N' or b == '+':
            continue
        length += 1
        if b == c:
            score += 1

    if test:
        print('search:')
        print(score, '/', length)
        print(barcode)
        print(code)

    return (length - score)

def find_matching_code (barcode, CODES, max_dist):
    code2dist = {}
    for k,v in CODES.items():
        if test:
            print(k)
        dist = barcode_dist(barcode, v)
        code2dist[k] = dist
    code2dist = dict(sorted(code2dist.items(), key=lambda item: item[1]))
    first_dist, second_dist, *_ = code2dist.values()
    if first_dist > second_dist: # smaller better, should not be equal
        print(code2dist)
        sys.exit('more than one CODE match with barcode')
    if first_dist <= max_dist:
        match_key =  next(iter(code2dist.keys()))
        if test:
            print(code2dist)
            print(match_key)
        return  match_key
    else:
        return False

def split_pe_fastq_by_barcode(input_fastq1, input_fastq2, output_directory):
    barcode_to_file = {}
    n_match = 0
    n_unknown = 0

    with gzip.open(input_fastq1, 'rt') as f1, gzip.open(input_fastq2, 'rt') as f2:
        while True:
            # Read four lines at a time for each read from both files
            header1 = f1.readline().strip()
            header2 = f2.readline().strip()
            if not header1 or not header2:
                break  # End of file
            
            sequence1 = f1.readline().strip()
            sequence2 = f2.readline().strip()
            plus_line1 = f1.readline().strip()
            plus_line2 = f2.readline().strip()
            quality1 = f1.readline().strip()
            quality2 = f2.readline().strip()
    
            if test:
                print(header1)
                print(header2)
                print(sequence1)
                print(sequence2)

            barcode1 = header1[-17:]
            barcode2 = header2[-17:]
            if barcode1 != barcode2:
                print(barcode1)
                print(barcode2)
                sys.exit('barcode1 NE barcode2') 
            else:
                barcode=barcode1

            key = find_matching_code(barcode, CODES, max_dist)

            if key: # match
                n_match += 1
                barcode_file1 = open(f"{output_directory}/{key}.R1.fastq", 'a')
                barcode_file2 = open(f"{output_directory}/{key}.R2.fastq", 'a')
            else: # not recognized
                n_unknown += 1
                barcode_file1 = open(f"{output_directory}/unknown/{barcode}.R1.fastq", 'a')
                barcode_file2 = open(f"{output_directory}/unknown/{barcode}.R2.fastq", 'a')
                if test:
                    print(barcode)
                    sys.exit('barcode not recognized')

            barcode_to_file[barcode] = (barcode_file1, barcode_file2)

            # Write the read pairs to the corresponding barcode files
            barcode_file1.write(f"{header1}\n{sequence1}\n{plus_line1}\n{quality1}\n")
            barcode_file2.write(f"{header2}\n{sequence2}\n{plus_line2}\n{quality2}\n")

    # Close all barcode-specific files
    print('n_match:', n_match)
    print('n_unkown:', n_unkown)
    print('n_total:', n_match + n_unkown)
    for barcode_file1, barcode_file2 in barcode_to_file.values():
        barcode_file1.close()
        barcode_file2.close()

if __name__ == "__main__":
    input_fastq1 = sys.argv[1]  # Replace with your input forward read FASTQ file
    input_fastq2 = sys.argv[2] # Replace with your input reverse read FASTQ file
    max_dist = int(sys.argv[3]) # smaller or equal, does not count N 
    output_directory = sys.argv[4]     # Replace with your desired output directory
    print('usage: python split.py test.R1.fastq.gz test.R2.fastq.gz 1 split_maxdist1')
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    if not os.path.exists(output_directory +'/'+ 'unknown'):
        os.mkdir(output_directory +'/'+ 'unknown')
    split_pe_fastq_by_barcode(input_fastq1, input_fastq2, output_directory)

