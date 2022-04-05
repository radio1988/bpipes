import os, re
files = os.listdir('fastq')
files = [file for file in files if re.search(r'R1.fastq.gz', file)]
SAMPLES = [file.replace('.R1.fastq.gz', '') for file in files]
for s in SAMPLES:
    print(s)
