blast7_file = 'gli.test.blast7'
min_len = 8

with open(blast7_file, 'r') as fp:
    line = fp.readline()
    i = 1
    while line:
        i += 1
        if line.startswith('#'):
            line = fp.readline()
            continue
        line = line.strip()
        [query_acc_ver, subject_acc_ver, identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score] = line.split('\t')
        effective_aligned_len = int(float(identity) /100 * int(alignment_length))
        bed = [subject_acc_ver, s_start, s_end, query_acc_ver, bit_score, '.', effective_aligned_len, -1, -1, -1]
        if effective_aligned_len >= min_len:
            bed = [str(x) for x in bed]
            print("\t".join(bed))
        line = fp.readline()
