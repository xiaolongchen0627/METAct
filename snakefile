configfile:"config/config.yaml"

import pandas as pd
#from snakemake.utils import validate
import os
#import pysam

units = pd.read_table(config["units"]).set_index("sample_name",drop=False)

import function.common_functions as cf

def is_paired(wildcards):
    if units.loc[wildcards.sample]["paired"] == "Yes" :
        return True
    else:
        return False

def Trimming(wildcards):
    if units.loc[wildcards.sample]["need_trim"] == "Yes" :
        return True
    else: 
        return False

def trim_inputs(wildcards):
    fq_input = [units.loc[wildcards.sample]["fq1"]]
    if(is_paired(wildcards)):
        fq_input.append(units.loc[wildcards.sample]["fq2"])
    return fq_input

def hg38_mapping_results(wildcards):
    hg38_input = [units.loc[wildcards.sample]["bam_2hg38"]]
    return hg38_input

rule all:  
    input:
        expand(os.path.join(config["working_dir"],"hg38","{sample}.hg38.bam.bai"),sample=units.sample_name)

rule trim_galore:
    input:
        trim_inputs
    output:
        os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz"),
        os.path.join(config["working_dir"],"trimmed","{sample}_val_2.fq.gz"),
        os.path.join(config["working_dir"],"trimmed","{sample}_trimmed.fq.gz")
    params:
        basename = "{sample}",
        trimmed_dir = config["working_dir"]+"/trimmed/",
        paired = is_paired
    run: 
        cmd = "trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} "
        if params.paired:
            cmd += "--paired {input} ; touch {params.trimmed_dir}/{params.basename}_trimmed.fq.gz"
        else :
            cmd += " {input}; touch {params.trimmed_dir}/{params.basename}_val_1.fq.gz;"
            cmd += "touch {params.trimmed_dir}/{params.basename}_val_2.fq.gz ;"
        shell(cmd)

rule align_2l1hs:
    input:
        fq1=os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz"),
        fq2=os.path.join(config["working_dir"],"trimmed","{sample}_val_2.fq.gz"),
        fqse=os.path.join(config["working_dir"],"trimmed","{sample}_trimmed.fq.gz")
    output:
        os.path.join(config["working_dir"],"l1hs","{sample}.sam")
    params:
        bt2_thread = config['bt2_thread'],
        bt2_idx = config['bt2_l1_idx'],
        paired = is_paired
    run:
        cmd = " bowtie2 --very-sensitive-local --thread {params.bt2_thread} -x {params.bt2_idx}"
        if params.paired:
            cmd +=" -U {input.fq1},{input.fq2} "
        else:
            cmd +=" -U {input.fqse} "
        cmd += "-S {output}"
        shell(cmd)

rule filter_l1hs_sam_editDistance:
    input:
        rules.align_2l1hs.output
    output:
        bam=os.path.join(config["working_dir"],"l1hs","{sample}.l1hs.bam"),
        index=os.path.join(config["working_dir"],"l1hs","{sample}.l1hs.bam.bai")
    params:
        editD = config['edit_distance']
    run:
        cmd = "awk '$1 ~ /^@/ || $0 ~ /NM:i:[0-{params.editD}]\t/ ' {input} | "
        cmd += "samtools view -Sb -F 4 - "
        cmd += " -o {output.bam} ;"
        cmd += "samtools sort {output.bam} -o {output.bam}.sort ; mv {output.bam}.sort {output.bam}; "
        cmd += "samtools index {output.bam}"
        shell(cmd)

rule extract_sclip2fq:
    input:
        bam=rules.filter_l1hs_sam_editDistance.output.bam
    output:
        fq=os.path.join(config["working_dir"],"l1hs","{sample}.sclip.fq")
    params:
        all_sclip = config['all_sclip'],
        sclip_length = config['min_sclip_len']
    run:
        cf.extractSoftclipReads(str(input.bam),str(output.fq),params.all_sclip,params.sclip_length)

rule align_2L1_sources:
    input :
        rules.extract_sclip2fq.output
    output: 
        fq=os.path.join(config["working_dir"],"l1hs","{sample}.None_L1.fq")
    params:
        bt2_idx = config['bt2_l1source_idx']
    run: 
        ### only keep the sclip reads that failed to map to other L1 sources
        cmd = " bowtie2 --very-sensitive-local -x {params.bt2_idx} -U {input} -S {output}.sam ;"
        cmd += " samtools view -Sb -f 4 {output}.sam -o {output}.none_l1.bam; "
        cmd += " picard SamToFastq -I {output}.none_l1.bam -F {output.fq} ;"
        shell(cmd)

rule realign_2hg38:
    input: 
        rules.align_2L1_sources.output
    output:
        os.path.join(config["working_dir"],"hg38","{sample}.hg38.sam")
    params:
        bt2_idx = config['bt2_hg38_idx']
    run:   
        cmd = " bowtie2 --very-sensitive-local -x {params.bt2_idx} -U {input} -S {output}"
        shell(cmd)

rule filter_hg38_sam:
    input:
        rules.realign_2hg38.output
    output:
        bam = os.path.join(config["working_dir"],"hg38","{sample}.hg38.bam"),
        bai = os.path.join(config["working_dir"],"hg38","{sample}.hg38.bam.bai")
    params:
        psl = config['L1_psl']
    run:
        cf.filter_hg38Sam(input,params.psl,str(output.bam)),
        cmd = 'samtools sort {output.bam} -o {output.bam}.sort ; mv {output.bam}.sort {output.bam}; '
        cmd += 'samtools index {output.bam} ' 
        shell(cmd)

rule score_sam:
    input:
        sclipSam = rules.filter_hg38_sam.output,
        refSam = hg38_mapping_results
    output:
        os.path.join(config["working_dir"],"hg38","{sample}.hg38.score.txt")
    run:
        cf.score_hg38Sam(sclipSam,refSam,str(output))
        


        

