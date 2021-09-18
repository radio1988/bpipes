import pandas as pd
from functools import reduce

class parse_meta_contrast_class:
  def __init__(self, y0, y1, y2, y3):
   self.contrast2contrast_name = y0
   self.contrast2treatmentSamples = y1
   self.contrast2controlSamples = y2
   self.contrasts = y3

def parse_meta_contrast (fmeta="meta.csv", fcontrast="contrast.csv"):
    """
    designed to get output dictionaries

    usage:
    parse_meta_contrast_obj=parse_meta_contrast(fmeta="meta.csv", fcontrast="contrast.csv") 
    output:
    # print("parse_meta_contrast_obj:", vars(parse_meta_contrast_obj))
    # {'contrast2contrast_name': {'contrast1': 'G1_vs_ctrl', 'contrast2': 'G2_vs_ctrl', 'contrast3': 'G1_G2_vs_ctrl'}, 
    # 'contrast2treatmentSamples': {'contrast1': ['2-1_S2', '2-2_S3', '2-3_S4'], 'contrast2': ['3-1_S5', '3-2_S6', '3-3_S7'], 'contrast3': ['2-1_S2', '2-2_S3', '2-3_S4', '3-1_S5', '3-2_S6', '3-3_S7']}, 
    # 'contrast2controlSamples': {'contrast1': ['1-2_S1'], 'contrast2': ['1-2_S1'], 'contrast3': ['1-2_S1']}}
    """
    ### parse meta
    meta = pd.read_csv(fmeta, comment='#')
    # todo: remove empty row/columns from meta.csv and contrast.csv
    #    sample group
    # 0  1-2_S1  ctrl
    # 1  2-1_S2    G1
    # 2  2-2_S3    G1
    # 3  2-3_S4    G1
    # 4  3-1_S5    G2
    # 5  3-2_S6    G2
    # 6  3-3_S7    G2

    ### group2sample
    temp = meta.copy()
    temp = temp[["sample", "group"]]
    temp['sample'] =  temp.groupby(['group']).transform(lambda x: ','.join(x))
    #                  sample group
    # 0                1-2_S1  ctrl
    # 1  2-1_S2,2-2_S3,2-3_S4    G1
    # 2  2-1_S2,2-2_S3,2-3_S4    G1
    # 3  2-1_S2,2-2_S3,2-3_S4    G1
    # 4  3-1_S5,3-2_S6,3-3_S7    G2
    # 5  3-1_S5,3-2_S6,3-3_S7    G2
    # 6  3-1_S5,3-2_S6,3-3_S7    G2
    temp['sample'] = temp['sample'].apply(lambda x: x.split(','))
    temp.index = temp['group']
    #                              sample group
    # group
    # ctrl                   [1-2_S1]  ctrl
    # G1     [2-1_S2, 2-2_S3, 2-3_S4]    G1
    # G1     [2-1_S2, 2-2_S3, 2-3_S4]    G1
    # G1     [2-1_S2, 2-2_S3, 2-3_S4]    G1
    # G2     [3-1_S5, 3-2_S6, 3-3_S7]    G2
    # G2     [3-1_S5, 3-2_S6, 3-3_S7]    G2
    # G2     [3-1_S5, 3-2_S6, 3-3_S7]    G2
    temp = temp.to_dict()
    # {'sample': {'ctrl': ['1-2_S1'], 'G1': ['2-1_S2', '2-2_S3', '2-3_S4'], 'G2': ['3-1_S5', '3-2_S6', '3-3_S7']}, 
    # 'group': {'ctrl': 'ctrl', 'G1': 'G1', 'G2': 'G2'}}
    group2sample=temp['sample']
    # {'ctrl': ['1-2_S1'], 'G1': ['2-1_S2', '2-2_S3', '2-3_S4'], 'G2': ['3-1_S5', '3-2_S6', '3-3_S7']}

    ### parse contrast
    contrast = pd.read_csv(fcontrast, comment='#')
    #         type treatment control
    # 0  contrast1        G1    ctrl
    # 1  contrast2        G2    ctrl
    # 2  contrast3     G1;G2    ctrl
    ### contrast2groups
    contrasts =  contrast['type'].tolist()  
    # ['contrast1', 'contrast2', 'contrast3']
    contrastNames = contrast.loc[:,['treatment', 'control']].agg("_vs_".join, axis=1).tolist()
    contrastNames = [x.replace(';', '_') for x in contrastNames] 
    # ['G1_vs_ctrl', 'G2_vs_ctrl', 'G1_G2_vs_ctrl']
    treatmentGroups = [x.split(';') for x in contrast['treatment']]
    # [['G1'], ['G2'], ['G1', 'G2']]
    controlGroups = [x.split(';') for x in contrast['control']]
    # [['ctrl'], ['ctrl'], ['ctrl']
    contrast2contrast_name=dict(zip(contrasts, contrastNames))
    # {'contrast1': 'G1_vs_ctrl', 'contrast2': 'G2_vs_ctrl', 'contrast3': 'G1_G2_vs_ctrl'}
    contrast2treatmentGroups=dict(zip(contrasts, treatmentGroups))
    # {'contrast1': ['G1'], 'contrast2': ['G2'], 'contrast3': ['G1', 'G2']}
    contrast2controlGroups=dict(zip(contrasts, controlGroups))
    # {'contrast1': ['ctrl'], 'contrast2': ['ctrl'], 'contrast3': ['ctrl']}

    ### contrast2groups2samples
    # input: contrast2treatmentGroups, contrast2controlGroups, group2sample
    # output: contrast2treatmentSamples, contrast2controlSamples
    contrast2treatmentSamples = {}
    for c, gs in contrast2treatmentGroups.items():
        # print(c, gs)
        # contrast3 ['G1', 'G2']
        ss = list(map(group2sample.get, gs))
        # [['2-1_S2', '2-2_S3', '2-3_S4'], ['3-1_S5', '3-2_S6', '3-3_S7']]
        ss = reduce(lambda x,y: x+y, ss)
        # ['2-1_S2', '2-2_S3', '2-3_S4', '3-1_S5', '3-2_S6', '3-3_S7']
        contrast2treatmentSamples[c] = ss
        # {'contrast1': ['2-1_S2', '2-2_S3', '2-3_S4'], 
        # 'contrast2': ['3-1_S5', '3-2_S6', '3-3_S7'], 
        # 'contrast3': ['2-1_S2', '2-2_S3', '2-3_S4', '3-1_S5', '3-2_S6', '3-3_S7']}
    contrast2controlSamples = {}
    for c, gs in contrast2controlGroups.items():
        ss = list(map(group2sample.get, gs))
        ss = reduce(lambda x,y: x+y, ss)
        contrast2controlSamples[c] = ss


    return parse_meta_contrast_class(contrast2contrast_name, 
                contrast2treatmentSamples, contrast2controlSamples, contrasts) # dictionaries in class

def get_treatment_bams_from_contrast(contrast="contrast1", root="results/clean_reads/", o = "parse_meta_contrast_obj"):
    treatment_samples = o.contrast2treatmentSamples[contrast]
    treatment_bams = list(map(lambda x:  root + str(x) + ".bam", treatment_samples))
    # print ("treatment_bams now:", treatment_bams)
    # treatment_bams now: ['clean_reads/2-1_S2.bam', 'clean_reads/2-2_S3.bam', 'clean_reads/2-3_S4.bam', 'clean_reads/3-1_S5.bam', 'clean_reads/3-2_S6.bam', 'clean_reads/3-3_S7.bam']
    #treatment_bams = reduce(lambda x, y: x+" "+y , treatment_bams)
    return treatment_bams

def get_control_bams_from_contrast(contrast="contrast1", root="results/clean_reads/", o = "parse_meta_contrast_obj"):
    control_samples = o.contrast2controlSamples[contrast]
    control_bams = list(map(lambda x:  root + str(x) + ".bam", control_samples))
    return control_bams

def get_contrast_name_from_contrast(contrast="contrast1", o = "parse_meta_contrast_obj"):
    # print("name", contrast, o.contrast2contrast_name[contrast])
    return o.contrast2contrast_name[contrast]

def get_narrowPeak_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    treat_pileup_bdg_names = []
    for contrast,name in zip(contrasts, contrast_names):
        treat_pileup_bdg_names.append("results/narrow_peaks_contrast_level/"+contrast+"/"+name+"_clean.narrowPeak")
    return treat_pileup_bdg_names

def get_broadPeak_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    treat_pileup_bdg_names = []
    for contrast,name in zip(contrasts, contrast_names):
        treat_pileup_bdg_names.append("results/broad_peaks_contrast_level/"+contrast+"/"+name+"_clean.broadPeak")
    return treat_pileup_bdg_names

def get_treat_pileup_bw_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [narrow_peaks_contrast_level/contrast1/G1_vs_ctrl_treat_pileup.bdg, 
    narrow_peaks_contrast_level/contrast2/G2_vs_ctrl_treat_pileup.bdg, 
    narrow_peaks_contrast_level/contrast3/G1_G2_vs_ctrl_treat_pileup.bdg]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    treat_pileup_bdg_names = []
    for contrast,name in zip(contrasts, contrast_names):
        treat_pileup_bdg_names.append(
            "results/narrow_peaks_contrast_level/"+contrast+"/"+name+"_treat_pileup.bw")
    return treat_pileup_bdg_names

def get_control_lambda_bdg_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [narrow_peaks_contrast_level/contrast1/G1_vs_ctrl_control_lambda.bdg, 
    narrow_peaks_contrast_level/contrast2/G2_vs_ctrl_control_lambda.bdg, 
    narrow_peaks_contrast_level/contrast3/G1_G2_vs_ctrl_control_lambda.bdg]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    control_lambda_bdg_names = []
    for contrast,name in zip(contrasts, contrast_names):
        control_lambda_bdg_names.append(
            "results/narrow_peaks_contrast_level/"+contrast+"/"+name+"_control_lambda.bdg")
    return control_lambda_bdg_names

def get_control_lambda_bw_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [narrow_peaks_contrast_level/contrast1/G1_vs_ctrl_control_lambda.bdg, 
    narrow_peaks_contrast_level/contrast2/G2_vs_ctrl_control_lambda.bdg, 
    narrow_peaks_contrast_level/contrast3/G1_G2_vs_ctrl_control_lambda.bdg]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    control_lambda_bdg_names = []
    for contrast,name in zip(contrasts, contrast_names):
        control_lambda_bdg_names.append("results/narrow_peaks_contrast_level/"+contrast+"/"+name+"_control_lambda.bw")
    return control_lambda_bdg_names


def get_narrow_count_names_from_contrasts(contrasts=["contrast1", "contrast2"], o = "parse_meta_contrast_obj"):
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append("results/narrow_peaks_contrast_level/"+contrast+"/"+name+"_count.txt")
    return outnames


def get_broad_count_names_from_contrasts(contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append("results/broad_peaks_contrast_level/"+contrast+"/"+name+"_count.txt")
    return outnames


def get_meme_peak_outname_from_contrasts(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [narrow_peaks_contrast_level/contrast1/memeSummit.2000/G1_vs_ctrl.finished
    ...
    ]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            "results/narrow_peaks_contrast_level/"+contrast\
            +"/meme_clean.real_peaks/"+name+'.finished')
        outnames.append(
            "results/broad_peaks_contrast_level/"+contrast\
            +"/meme_clean.real_peaks/"+name+'.finished')
    return outnames



def get_meme_summit_outname_from_contrasts(
    contrasts=["contrast1", "contrast2"], 
    PEAK_WIDTH=["100","1000", "2000"], 
    o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [narrow_peaks_contrast_level/contrast1/memeSummit.2000/G1_vs_ctrl.finished
    ...
    ]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        for w in PEAK_WIDTH:
            outnames.append(
                "results/narrow_peaks_contrast_level/"+contrast\
                +"/memeSummit." + str(w) + "/" + name + '.finished') # broadpeak no summit
    return outnames


def get_meme_summit_split_outname_from_contrasts(
    contrasts=["contrast1", "contrast2"], 
    PEAK_WIDTH=["100","1000", "2000"], 
    CHRS=['chr1', 'chrX', 'chrY'],
    o = "parse_meta_contrast_obj"):
    """
    Learn: Good trick to use tagets input to do contrast2contrast_name and more

    output example:
    [
    narrow_peaks_contrast_level/contrast1/memeSummit_chr.100_chrX/G1_vs_ctrl.finished
    ...
    ]
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        for w in PEAK_WIDTH:
            for chr in CHRS:
                outnames.append(
                    "results/narrow_peaks_contrast_level/"+contrast\
                    +"/memeSummit_chr." + str(w) + "_" + chr + '/' + name + '.finished')
    return outnames


def get_signalHeatmap_ContrastPeak_outname(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    output example:
    "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak.pdf"
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            'results/narrow_peaks_contrast_level/'+contrast+'/'+name+'_clean.real.narrowPeak.pdf')
        outnames.append(
            'results/broad_peaks_contrast_level/'+contrast+'/'+name+'_clean.real.broadPeak.pdf')
    return outnames

def get_corrHeatmap_peakCount_outname(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    output example:
    "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.heatmap.pdf"
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            'results/narrow_peaks_contrast_level/'+contrast+'/'+name+'_count.heatmap.pdf')
        outnames.append(
            'results/broad_peaks_contrast_level/'+contrast+'/'+name+'_count.heatmap.pdf')
    return outnames


def get_chippeakanno_outname(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    output example:
    "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx"
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            'results/narrow_peaks_contrast_level/'+contrast\
            +'/'+name+'_clean.real.narrowPeak.final_anno.xlsx')
        outnames.append(
            'results/broad_peaks_contrast_level/'+contrast\
            +'/'+name+'_clean.real.broadPeak.final_anno.xlsx')
    return outnames



def get_cpm_filter_outname(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    output example:
        PEAK_UP="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak",,
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            'results/narrow_peaks_contrast_level/'+contrast+'/'+name+'_clean.real.narrowPeak')
        outnames.append(
            'results/broad_peaks_contrast_level/'+contrast+'/'+name+'_clean.real.broadPeak')
    return outnames


def get_enrichment_analysis_outname(
    contrasts=["contrast1", "contrast2"], 
    o = "parse_meta_contrast_obj"):
    """
    output example:
    "results/{narrowbroad}_peaks_contrast_level/" \
                 + "{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak.enrichment.finished"
    """
    contrast_names = map(o.contrast2contrast_name.get, contrasts)
    outnames = []
    for contrast,name in zip(contrasts, contrast_names):
        outnames.append(
            'results/narrow_peaks_contrast_level/'+contrast+\
            '/'+name+'_clean.real.narrowPeak.enrichment.finished')
        outnames.append(
            'results/broad_peaks_contrast_level/'+contrast+\
            '/'+name+'_clean.real.broadPeak.enrichment.finished')
    return outnames

    
