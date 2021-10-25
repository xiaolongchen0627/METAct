configfile:"config/config.yaml"

import pandas as pd
from snakemake.utils import validate
import os
#import pysam

units = pd.read_table(config["units"]).set_index("sample_name",drop=False)

from function.common_functions import extractSoftclipReads

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

rule all:  
    input:
        expand(os.path.join(config["working_dir"],"hg38","{sample}.hg38.sam"),sample=units.sample_name)

rule trim_galore:
    input:
        trim_inputs
    output:
        temp(os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz")),
        temp(os.path.join(config["working_dir"],"trimmed","{sample}_val_2.fq.gz")),
        temp(os.path.join(config["working_dir"],"trimmed","{sample}_trimmed.fq.gz"))
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
        cmd += "samtools sort {output.bam} {output.bam}.sort ; mv {output.bam}.sort.bam {output.bam}; "
        cmd += "samtools index {output.bam}"
        print(cmd)
        shell(cmd)

rule extract_sclip2fq:
    input:
        bam=rules.filter_l1hs_sam_editDistance.output.bam
    output:
        fq=temp(os.path.join(config["working_dir"],"l1hs","{sample}.sclip.fq"))
    run:
        extractSoftclipReads(str(input.bam),str(output.fq))

rule align_2L1_sources:
    input :
        rules.extract_sclip2fq.output
    output: 
        fq=temp(os.path.join(config["working_dir"],"l1hs","{sample}.None_L1.fq"))
    params:
        bt2_idx = config['bt2_l1source_idx']
    run: 
        ### only keep the sclip reads that failed to map to other L1 sources
        cmd = " bowtie2 --very-sensitive-local -x {params.bt2_idx} -U {input} -S {output}.sam ;"
        cmd += " samtools view -Sb -f 4 {output}.sam -o {output}.none_l1.bam; "
        cmd += " picard SamToFastq -I {output}.none_l1.bam -F {output.fq} ;"
        cmd += " rm {output}.none_l1.bam {output}.sam ;"
        shell(cmd)

rule align_2hg38:
    input: 
        rules.align_2L1_sources.output
    output:
        os.path.join(config["working_dir"],"hg38","{sample}.hg38.sam")
    params:
        bt2_idx = config['bt2_hg38_idx']
    run:   
        cmd = " bowtie2 --very-sensitive-local -x {params.bt2_idx} -U {input} -S {output}"
        shell(cmd)

#rule filter_hg38_sam_editDistance:
#    input:
#        rules.align_2hg38.output:
#    output:
#        os.path.join(config["working_dir"],"hg38","{sample}.hg38.bam")
#    params:
#        psl = config['L1_psl']
#    run:
#        cmd = " awk '($1 ~ /_Ci:.{5,6}[SM]$/ && $0 ~ /XM:i:[0-1]/ )' {input} |"
#        cmd += " | awk 'NR==FNR&&NR>=6
        cmd += "{split($10,name,".");"
        cmd += "chr[$10]=name[1];s[$10]=name[2];e[$10]=name[3];n[$10]=name[4];ls[$10]=$15;lend[$10]=$16;print(le); next}"
        # psl file , blat results of L1PA1 to L1PA7 towards L1HS 
               chr[$10]=name[1];s[$10]=name[2];e[$10]=name[3];n[$10]=name[4];ls[$10]=$15;lend[$10]=$16;print(le); next} "


#awk '($1 ~ /_Ci:.{5,6}[SM]$/ && $0 ~ /XM:i:[0-1]/ )' SJBT031173_D1.hg38.sam | awk 'NR==FNR&&NR>=6{split($10,name,".");chr[$10]=name[1];s[$10]=name[2];e[$10]=name[3];n[$10]=name[4];ls[$10]=$15;lend[$10]=$16;print(le); next} {split($1,tmp,"pos:");split(tmp[2],pos,"_Ci");for(each in chr){if(chr[each]==$3&&$4>=s[each]-10000&&$4<=e[each]+10000&&pos[1]>=ls[each]-200&&pos[1]<=lend[each]+200)print $0"\t"each}}' ~/FOXR2/rmsk.L1HS.consensus.psl -  |grep -v "^$" > tmp1

