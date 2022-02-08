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

def need_trim(wildcards):
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
        expand(os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz"),sample=units.sample_name)

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
        paired = is_paired,
        trimmed = need_trim
    run: 
        cmd = ""
        if params.trimmed:
            cmd += "trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} "
            if params.paired:
                cmd += "--paired {input} ; touch {params.trimmed_dir}/{params.basename}_trimmed.fq.gz"
            else :
                cmd += " {input}; touch {params.trimmed_dir}/{params.basename}_val_1.fq.gz;"
                cmd += "touch {params.trimmed_dir}/{params.basename}_val_2.fq.gz ;"
        else:
            
            if params.paired:
                if input[0].endswith('.gz'):                                        
                    cmd += " ln -s {input[0]} {params.trimmed_dir}/{params.basename}_val_1.fq.gz; "
                    cmd += " ln -s {input[1]} {params.trimmed_dir}/{params.basename}_val_2.fq.gz; "
                else:                    
                    cmd += " gzip -c {input[0]} > {params.trimmed_dir}/{params.basename}_val_1.fq.gz;"
                    cmd += " gzip -c {input[1]} > {params.trimmed_dir}/{params.basename}_val_2.fq.gz;"
                cmd += " touch {params.trimmed_dir}/{params.basename}_trimmed.fq.gz; "
            else:
                if input[0].endswith('.gz'):
                    cmd += " ln -s {input[0]} {params.trimmed_dir}/{params.basename}_trimmed.fq.gz;"
                else : 
                    cmd += " gzip -c {input[0]} > {params.trimmed_dir}/{params.basename}_trimmed.fq;"
                cmd += "touch {params.trimmed_dir}/{params.basename}_val_1.fq.gz;"
                cmd += "touch {params.trimmed_dir}/{params.basename}_val_2.fq.gz;"

        shell(cmd )

rule align2CustomGenome:
    input:
        fq1=os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz"),
        fq2=os.path.join(config["working_dir"],"trimmed","{sample}_val_2.fq.gz"),
        fqse=os.path.join(config["working_dir"],"trimmed","{sample}_trimmed.fq.gz")
    output:
        os.path.join(config["working_dir"],config["customGenomeName"],"{sample}.sam")
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

rule filter_sam_editDistance:
    input:
        rules.align2CustomGenome.output
    output:
        bam=os.path.join(config["working_dir"],config["customGenomeName"],"{sample}.",config["customGenomeName"],".bam"),
        bai=os.path.join(config["working_dir"],config["customGenomeName"],"{sample}.",config["customGenomeName"],".bam.bai")
    params:
        editD = config['edit_distance']
    run:
        cmd = "awk '$1 ~ /^@/ || $0 ~ /NM:i:[0-{params.editD}]\t/ ' {input} | "
        cmd += "samtools view -Sb -F 4 -q 20 - "
        cmd += " -o {output.bam} ;"
        cmd += "samtools sort {output.bam} -o {output.bam}.sort ; mv {output.bam}.sort {output.bam}; "
        cmd += "samtools index {output.bam}"
        shell(cmd)

rule extract_sclip2fq:
    input:
        bam=rules.filter_sam_editDistance.output.bam
    output:
        fq=os.path.join(config["working_dir"],config["customGenomeName"],"{sample}.sclip.fq")
    params:
        sclip_length = config['min_sclip_len']
    run:
        cf.extractSoftclipReads(str(input.bam),str(output.fq),params.sclip_length)


rule realign_2hg38:
    input: 
        rules.extract_sclip2fq.output
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
        psl = config['L1_psl'],
        focus_site = config['interestingSite'],
        target = config['customGenomeName']
    run:
        cf.filter_hg38Sam(input,params.psl,str(output.bam),params.focus_site,params.target),
        cmd = 'samtools sort {output.bam} -o {output.bam}.sort ; mv {output.bam}.sort {output.bam}; '
        cmd += 'samtools index {output.bam} ' 
        shell(cmd)

rule score_sam:
    input:
        sclipSam = rules.filter_hg38_sam.output.bam,
        refSam = hg38_mapping_results
    output:
        sam = os.path.join(config["working_dir"],"hg38","{sample}.hg38.rev.sam"),
        fqname = os.path.join(config["working_dir"],"hg38","{sample}.hg38.rev.fq.name.txt")
    run:
        cf.scoring_alignment(input.refSam,input.sclipSam,str(output.sam))
        
rule get_fastq_for_validation:
    input :
        fq1=rules.align2CustomGenome.input.fq1,
        fq2=rules.align2CustomGenome.input.fq2,
        fqname=rules.score_sam.output.fqname
    params:
        paired = is_paired
    output: 
        fq1 = os.path.join(config["working_dir"],"hg38","{sample}.hg38.rev.r1.fq"),
        fq2 = os.path.join(config["working_dir"],"hg38","{sample}.hg38.rev.r2.fq")
    run:
        cmd = ' '
        if params.paired:
            cmd +=" zcat {input.fq1} | grep -f {input.fqname} -A 3 -h  --no-group-separator > {output.fq1};"
            cmd += " zcat  {input.fq2}  | grep -f {input.fqname} -A 3 -h  --no-group-separator> {output.fq2};"
        else:
            cmd += "touch {output.fq1} ;"
            cmd += "touch {output.fq2} "
        shell(cmd)


        

