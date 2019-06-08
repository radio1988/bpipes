import os
from os import path
from itertools import compress  # slice list with boolean idx
import re
import random
import json
import subprocess


# parse args
import argparse

parser = argparse.ArgumentParser(
    description='''
    usage1: fix_fusion_transcripts.py -query cuffmerge.gtf -db1 refseq.gtf -db2 ensembl.gtf -mode find --min-overlap 0.0
    usage2: fix_fusion_transcripts.py -query cuffmerge.gtf -db1 refseq.gtf -db2 ensembl.gtf -mode fix (after usage1)
    1. find transcripts in query.gtf that overlaps db1.gtf or db2.gtf, and create query.fusion.black.tr_list
    2. find transcripts in db1.gtf that overlaps black-query-transcript, and output db1.tba.tr_list
    3. if not 2, find transcripts in db2.gtf that overlaps, and output db2.tba.tr_list
    ''')
parser.add_argument('-query', help='query gtf, e.g. cuffmerge.gtf')
parser.add_argument('-db1',  help='database gtf with higher priority, e.g. refseq.gtf')
parser.add_argument('-db2',  help='database gtf with lower priority, e.g. ensembl.gtf')
parser.add_argument('-mode', default=None, help='find: find from db and create json; fix: read json and fix gtf')
parser.add_argument('--min-overlap', default=0.0, type=float, help="min overlap required to call fusion in find mode")
parser.add_argument('--chr-length', default=None,
                    help='if specified, fix lines in gtf that exceeds the boundries of chrs'
                    'e.g. GRZ11.ucsc.chr_length.txt: chr length')

args = parser.parse_args()

# check args
if(path.exists(args.query) & path.exists(args.db1) & path.exists(args.db2)):
    print (args)
    print ('All input files exist on disk')
else:
    raise Exception('Error: Some of', args, "Don't exist on disk")

if args.mode is None:
    raise Exception ("Error: Please specify a mode")


# functions

# GTF processing
def parse_gtf_line(line):
    """
    None for header and lower case type
    :param line:
    :return: gtf_dict (chr, name, type, start, end, anno, chrstrand)
    """
    line = line.rstrip()
    if re.match(r'^#', line):
        return None
    (chr, name, type, start, end, placer, strand, placer2, anno) = line.split('\t')
    anno = anno.replace(" ", "")
    chrstrand = chr + strand

    type = type.lower()

    gtf_dict = {"chr":chr,
                "name":name,
                "type":type,
                "start":start,
                "end":end,
                "anno":anno,
                "chrstrand":chrstrand}

    return gtf_dict


def read_gtf(fname):
    """
    skip header and type=gene lines
    :param fname: gtf
    :return: tr2exons, tr2gene (dictionaries)
    """
    tr2gene = {}
    tr2exons = {}
    tr2chrstrand = {}
    chrstrand2tr = {}
    gene2tr = {}
    print("reading ", fname)
    with open(fname, "r") as f:
        for line in f:
            line = line.rstrip()
            # print(line)

            if re.match(r'^#', line):
                continue

            (chr, name, type, start, end, placer, strand, placer2, anno) = line.split('\t')
            anno = anno.replace(" ", "")
            chrstrand = chr + strand

            type=type.lower()
            if type == "gene":
                continue

            p = re.compile("transcript_id([^;]*);")
            transcript_id = p.findall(anno)[0].replace("\"", "")

            p = re.compile("gene_id([^;]*);")
            p1 = re.compile("gene_id")
            if re.search("gene_id", anno) is None:
                continue  # skip lines without gene_id for refseq.gtf
            gene_id = p.findall(anno)[0].replace("\"", "")

            exon = [int(start), int(end)]
            # print(chr, strand, exon, transcript_id, gene_id)
            tr2gene.update({transcript_id: gene_id})
            gene2tr.update({gene_id: transcript_id})
            tr2chrstrand.update({transcript_id: chr+strand})

            if transcript_id in tr2exons.keys():
                tr2exons[transcript_id].append(exon)
            else:
                tr2exons.update({transcript_id: [exon]})

            if chrstrand in chrstrand2tr.keys():
                if chrstrand2tr[chrstrand] != transcript_id:
                    chrstrand2tr[chrstrand].append(transcript_id)
            else:
                chrstrand2tr.update({chrstrand: [transcript_id]})

        with open(fname+'.g2t', 'w') as outfile:
            for i, gene in enumerate(gene2tr):
                outline = str(gene) + "," + str(gene2tr[gene]) + "\n"
                outfile.write(outline)

        with open(fname + '.t2g', 'w') as outfile:
            for i, tr in enumerate(tr2gene):
                outline = str(tr) + "," + str(tr2gene[tr]) + "\n"
                outfile.write(outline)



    return tr2exons, tr2gene, tr2chrstrand, chrstrand2tr

def test_gtf_reader(fname):
    tr2exons, tr2gene, tr2chrstrand, chrstrand2tr = read_gtf(fname)
    # tr_key = random.choice(1) in tr2exons.keys()
    # chr_key = chrstrand2tr.keys()[0]
    # print(chrstrand2tr[chr_key])
    # print(tr2exons[tr_key])

    for key, value in tr2exons.items():
        print(key, value)
        break

    for key, value in tr2gene.items():
        print(key, value)
        break

    for key, value in tr2chrstrand.items():
        print(key, value)
        break

    for key, value in chrstrand2tr.items():
        print(key, value)
        break

def gtf_line_match_string_list(line, strlist):
    fusion_tr_match_list = [re.search(x, line) for x in strlist]
    no_match = all(x is None for x in fusion_tr_match_list)
    if no_match:
        return False
    else:
        matched_gene_list = list(compress(strlist, fusion_tr_match_list))
        return matched_gene_list

def get_transcript_id(line):
    if re.match(r'^#', line):
        return None
    line = line.replace(" ", "")
    p = re.compile("transcript_id([^;]*);")
    if p.search(line) is None:
        return None  # skip lines without gene_id for refseq.gtf

    transcript_id = p.findall(line)[0].replace("\"", "")
    return transcript_id

def get_gene_id(line):
    if re.match(r'^#', line):
        return None

    line = line.replace(" ", "")
    p = re.compile("gene_id([^;]*);")
    if p.search(line) is None:
        return None  # skip lines without gene_id for refseq.gtf

    gene_id = p.findall(line)[0].replace("\"", "")
    return gene_id

def append_write_gtf(line, fname):
    with open(fname, "a") as outfile:
        outfile.write(line)

def sort_exons(exons):
    """
    sort_exons(tr2exons['tr1'])
    the input object exons will be altered
    """
    exons.sort(key=lambda x: x[0])

def exon_length(exons):
    """
    total exon length of transcript
    """
    sort_exons(exons)
    return sum([x[1]-x[0]+1 for x in exons])  # gtf inclusive ranges => +1

def transcript_range(exons):
    """
    sorting included
    """
    sort_exons(exons)
    return [exons[0][0], exons[-1][1]]


def transcript_length(exons):
    """
    genomic length of the transcript
    """
    range_=transcript_range(exons)
    return range_[1] - range_[0] + 1


def exon_overlap_bp(exon1, exon2):
    """
    Only input two exon, not two exons
    :param exon1: e.g.[1001, 1100]
    :param exon2: e.g.[1011, 1110]
    :return: e.g. 90
    """
    return min(exon1[1], exon2[1]) - max(exon1[0], exon2[0]) + 1

def transcripts_overlap_fraction(exons1, exons2):
    """
    sorting included
    overlap_fraction = max(overlap_bp/length1, overlap_bp/length2)
    transcript length now
    :param exons1: for a transcript
    :param exons2:
    :return: overlap_fraction in float [0-1]
    """
    # quick screening (transcript level)
    transcript_range1 = transcript_range(exons1)
    transcript_range2 = transcript_range(exons2)
    if transcript_range1[1] < transcript_range2[0]:
        return 0
    if transcript_range1[0] > transcript_range2[1]:
        return 0

    # todo: accurate (exon level), now transcript level
    min_transcript_length = min(transcript_length(exons1), transcript_length(exons2))
    overlap_fraction = exon_overlap_bp(transcript_range1, transcript_range2) / min_transcript_length
    return round(overlap_fraction, 3)

def each_query_db_hits(query_exons, db_tr2exons, min_overlap_fraction):
    """
    :param query_exons:
    :param db_tr2exons:
    :return: db_transcript_ids that overlaps query_transcript
    """
    db_transcripts = db_tr2exons.keys()
    fractions = [transcripts_overlap_fraction(query_exons, db_tr2exons[x]) for x in db_transcripts]
    idx = [fraction > min_overlap_fraction for fraction in fractions]
    return list(compress(db_transcripts, idx))  # list of transcript ids ['tr1', 'tr2']

def queries_db_hits(query_tr2exons, db_tr2exons, min_overlap):
    """
    :param query_tr2exons:
    :param db_tr2exons:
    :param min_overlap:
    :return: [query.fusion.black.tr_list, db.add.white.tr_list]
    [query-black-0, []]  # no hit
    [query-black-1, [tr-white-1a, tr-white-1b]]
    [query-black-2, [tr-white-2a, tr-white-2b, tr-white-2c]]
    """
    query_transcripts = query_tr2exons.keys()
    db_hits = [each_query_db_hits(query_tr2exons[x], db_tr2exons, min_overlap) for x in query_transcripts]
    return list(zip(query_transcripts, db_hits))

def find_fusion_genes(full_res, db_tr2gene):
    quries = [x[0] for x in full_res]
    hits_genes = [[db_tr2gene[y] for y in x[1]] for x in full_res]
    hits_genes = [sorted(set(x)) for x in hits_genes]
    idx = [len(x) > 1 for x in hits_genes]  # more than one db_genes hit
    fusion_query_transcripts= list(compress(quries, idx))
    fution_db_genes = list(compress(hits_genes, idx))
    return list(zip(fusion_query_transcripts, fution_db_genes)) # [('tr1', ['Gene1', 'Gene2']), ('tr2', ['Gene1',
    # 'Gene2', 'Gene3'])

def write_fusion_res(fusion_res, label="query_db_hits"):
    with open(label+".csv", "a") as outfile:
        for key, value in fusion_res:
            out_txt = key + "," + ",".join(value) + "\n"
            outfile.write(out_txt)

    # with open(label+'.json', 'a') as outfile:
    #     json.dump(fusion_res, outfile)

def test_logic_code(min_overlap = 0.05):
    # build test data
    exons1 = [[1001, 1100], [100, 110], [1, 10], [30, 40]]
    exons2 = [[1011, 1110], [2011, 2110]]
    exons3 = [[2011, 2110], [3011, 3110], [4011, 4110]]
    exons4 = [[8011, 8110], [9011, 9110]]

    query_tr2gene_id = {'tr1':"gene1", 'tr2':'gene2', 'tr4':'gene4'}
    db_tr2gene = {'TR1':"Gene1", 'TR2':'Gene2', 'TR3':'Gene3'}

    query_tr2exons = {'tr1':exons1, 'tr2':exons2, 'tr4':exons4} # (tr2exons['tr1'][0][0]) # tr1, first exon, exon 0,
    # start-site
    db_tr2exons =  {'TR1':exons1, 'TR2':exons2, 'TR3':exons3}
    tr2strand = {'tr1':"+", 'tr2':'+'}
    query_transcripts = ['tr1', 'tr2']
    db_transcripts = ['tr1', 'tr2']

    print("Exons before soring:", exons1) # [[1001, 1100], [1, 100]]
    sort_exons(query_tr2exons['tr1'])
    print("Exons after soring:", exons1)  # [[1, 100], [1001, 1100]]

    print("Exon length for tr1:", exon_length(query_tr2exons['tr1']))  # 200

    print("transcript range:", transcript_range(query_tr2exons['tr1']))  # [1, 1100]

    print("exon_overlap_bp", exon_overlap_bp(exons1[-1], exons2[0]))  # 90
    print("transcripts_overlap_fraction", transcripts_overlap_fraction(exons1, exons2))  # 0.082

    # test finding overlap for all query transcripts
    # for one query
    db_transcript_hits = each_query_db_hits(exons1, db_tr2exons, min_overlap)
    print("db transcript hits for tr1:", db_transcript_hits)  #  ['TR1', 'TR2']

    # for all quries
    print ("query_tr2exons", query_tr2exons)
    print ("db_tr2exons", db_tr2exons)
    full_res = queries_db_hits(query_tr2exons, db_tr2exons, min_overlap)
    print("all query-db hits:", full_res)
    # [('tr4', []), ('tr1', ['TR2', 'TR1']), ('tr2', ['TR2', 'TR1', 'TR3'])]

    # find fusion genes
    fusion_res = find_fusion_genes(full_res, db_tr2gene)
    print("fusion hits:", fusion_res)

def find_mode():
    # read
    query_tr2exons, query_tr2gene, query_tr2chrstrand, query_chrstrand2tr = read_gtf(args.query)
    db1_tr2exons, db1_tr2gene, db1_tr2chrstrand, db1_chrstrand2tr = read_gtf(args.db1)
    db2_tr2exons, db2_tr2gene, db2_tr2chrstrand, db2_chrstrand2tr = read_gtf(args.db2)

    # find overlap
    # for each chrstrand
    # clear previous output
    try:
        [os.remove(x) for x in ["query_db1_hits.csv", "query_db2_hits.csv", "query_db1_hits.json", "query_db2_hits.json"]]
    except (FileNotFoundError):
        pass

    for idx, chrstrand in enumerate(sorted(query_chrstrand2tr)):
        # query
        query_transcript_list = query_chrstrand2tr[chrstrand]
        print(idx+1, '/', len(query_chrstrand2tr.keys()), 'chrstrand', chrstrand)

        query_tr2exons_ = {x:query_tr2exons[x] for x in query_transcript_list}
        print(len(query_tr2exons_.keys()), '/', len(query_tr2exons.keys()), 'query transcripts')

        # db1
        try:
            db1_transcript_list = db1_chrstrand2tr[chrstrand]
            db1_tr2exons_ = {x: db1_tr2exons[x] for x in db1_transcript_list}
            print(len(db1_tr2exons_.keys()), '/', len(db1_tr2exons.keys()), 'db1 transcripts')
            full_res1 = queries_db_hits(query_tr2exons_, db1_tr2exons_, args.min_overlap)
            fusion_res1 = find_fusion_genes(full_res1, db1_tr2gene)
            write_fusion_res(fusion_res1, label="query_db1_hits")
        except KeyError:
            print(args.db1, "does not have", chrstrand)

        # db2
        try:
            db2_transcript_list = db2_chrstrand2tr[chrstrand]
            db2_tr2exons_ = {x:db2_tr2exons[x] for x in db2_transcript_list}
            print(len(db2_tr2exons_.keys()), '/', len(db2_tr2exons.keys()), 'db2 transcripts')
            full_res2 = queries_db_hits(query_tr2exons_, db2_tr2exons_, args.min_overlap)
            fusion_res2 = find_fusion_genes(full_res2, db2_tr2gene)
            write_fusion_res(fusion_res2, label="query_db2_hits")
        except KeyError:
            print(args.db2, "does not have", chrstrand)

# fix mode functions
def read_saved_csv(fname):
    import csv
    query_transcripts = []
    db_genes = []
    # g2t = {}
    with open(fname) as file:
        readCSV = csv.reader(file, delimiter=',')
        for row in readCSV:
            tr = row[0]
            genes = row[1:None]
            query_transcripts.append(tr)
            db_genes.append(genes)
    return query_transcripts, db_genes


def read_chr_length(fname):
    chr2l = {}
    with open(fname, "r") as f:
        for line in f:
            line = line.rstrip()
            line = line.replace("\n", "")

            if re.match(r'^#', line):
                continue
            elif len(line) < 1:
                continue

            (chr, length) = line.split('\t')
            length = int(length)
            chr2l.update({chr: length})

    return chr2l


def fix_chr_length(gtf_fname, chr_length_fname):
    print("fix_chr_length for", gtf_fname, 'with', chr_length_fname)

    # read chr_length file
    if (path.exists(chr_length_fname)):
        print ("chr_length file:", chr_length_fname, 'exist on disk')
    else:
        raise Exception('Error: ', chr_length_fname, "Don't exist on disk")

    chr2l = read_chr_length(chr_length_fname)
    # [print(key, chr2l[key]) for key in chr2l.keys()]

    # read gtf
    print("reading ", gtf_fname)
    with open(gtf_fname + "fix", "w") as outf:
        with open(gtf_fname, "r") as f:
            for line in f:
                if re.match(r'^#', line):
                    outf.write(line)
                    continue
                (chr, name, type, start, end, placer, strand, placer2, anno) = line.split('\t')
                start=int(start)
                end=int(end)
                if start < 1:
                    start = 1
                if chr2l.get(chr) is None:
                    pass
                elif end >= chr2l[chr]:
                    print ("Line exceed chr length:", line)
                    end = chr2l[chr] - 1  # -1 just to be safe
                line = "\t".join([chr, name, type, str(start), str(end), placer, strand, placer2, anno])
                outf.write(line)

    rename_cmd = "mv  " + gtf_fname + "fix  " + gtf_fname
    print(rename_cmd)
    process = subprocess.Popen(rename_cmd, shell=True, executable='/bin/bash')

    return 1



def fix_mode():
    print("Fix mode:")

    # prep
    fusion_trs1, db1_genes = read_saved_csv("query_db1_hits.csv")
    fusion_trs2, db2_genes = read_saved_csv("query_db2_hits.csv")
    fusion_trs1 = set(fusion_trs1)
    fusion_trs2 = set(fusion_trs2)
    fusion_trs2_specific = fusion_trs2.difference(fusion_trs1)
    fusion_trs_both = fusion_trs1.union(fusion_trs2)
    print("fusion_trs1", len(fusion_trs1))
    print("fusion_trs2", len(fusion_trs2))
    print("fusion_trs2_specific", len(fusion_trs2_specific))
    print("fusion_trs_both", len(fusion_trs_both))
    db1_genes_set = set([item for sublist in db1_genes for item in sublist])
    db2_genes_set = set([item for sublist in db2_genes for item in sublist])
    db2_specific_genes_set = db2_genes_set.difference(db1_genes_set)

    # log
    with open('query_removed_transcripts.txt', 'w') as f:
        print("Dumping query_removed_transcripts.txt")
        query_removed_transcripts = list(fusion_trs_both)
        query_removed_transcripts.sort()
        json.dump(query_removed_transcripts, f)

    with open('db1_added_genes.txt', 'w') as f:
        print("Dumping db1_added_genes.txt")
        added_db1_genes = list(db1_genes_set)
        added_db1_genes.sort()
        json.dump(added_db1_genes, f)

    with open('db2_added_genes.txt', 'w') as f:
        print("Dumping db2_added_genes.txt")
        added_db2_genes = list(db2_specific_genes_set)
        added_db2_genes.sort()
        json.dump(added_db2_genes, f)

    # work
    # reset output
    try:
        gtf_files = ["query.kept.gtf", "query.removed.gtf", "db1.white.gtf", "db2.white.gtf"]
        [print("resetting", x) for x in gtf_files]
        [os.remove(x) for x in gtf_files]
    except (FileNotFoundError):
        pass

    # output
    with open(args.query) as f:
        print("keeping non-fusion transcripts (both db1 and db2) from query.gtf")
        for line in f:
            transcript_id = get_transcript_id(line)
            if transcript_id not in fusion_trs_both:
                append_write_gtf(line, 'query.kept.gtf')
            else:
                append_write_gtf(line, 'query.removed.gtf')

    with open(args.db1) as f:
        print('adding whitelist genes (all transcripts) from db1.gtf')
        for line in f:
            gene_id = get_gene_id(line)
            if gene_id in db1_genes_set:
                append_write_gtf(line, 'db1.white.gtf')

    with open(args.db2) as f:
        print('adding whitelist genes (all transcripts) from db2.gtf, that are not added in db1_whitelist_genes')
        for line in f:
            gtf_dict = parse_gtf_line(line)
            if gtf_dict is None:
                continue
            elif gtf_dict['type'] == 'gene':
                continue  # skip gene lines in ensembl
            else:
                gene_id = get_gene_id(line)
                if gene_id in db2_specific_genes_set:
                    append_write_gtf(line, 'db2.white.gtf')

    # sort
    for f in gtf_files:
        sort_gtf_cmd = "sort -k1,1 -k4,4n " + f + " > " + f + ".sorted.gtf " + \
                       "&& mv " + f + ".sorted.gtf " + f
        print(sort_gtf_cmd)
        process = subprocess.Popen(sort_gtf_cmd, shell=True, executable='/bin/bash')

    # fix exons with wrong ranges (introduced by cuffmerge's incompatibility with UCSC genome browser)
    if args.chr_length:
        [fix_chr_length(x, args.chr_length) for x in gtf_files]

    # merge, sort outputs
    merge_cmd = "cat query.kept.gtf db1.white.gtf db2.white.gtf |sort -k1,1 -k4,4n > fusion_fixed.gtf"
    print(merge_cmd)
    process = subprocess.Popen(merge_cmd, shell=True, executable='/bin/bash')







# Tests
# test_logic_code()
# test_gtf_reader(args.query)

# Work
if args.mode == 'find':
    find_mode()
elif args.mode == 'fix':
    fix_mode()
else:
    raise Exception("Error: unkown mode")

