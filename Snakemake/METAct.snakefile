"""
A pipeline to detect somatic mobile element driven activation transcription events in RNA-seq.
Author: Xiaolong Chen 10/13/2021
"""

import os


configfile: "config.yaml"
# configfile: "sample.config"

rule all:
	input: 
		r1 = "{wd}/{sample}/{sample}_val_1.fq.gz",
		r2 = "{wd}/{sample}/{sample}_val_2.fq.gz"


include: TrimGalore.snakefile

dependencies:






#rule bowtie2_map_step1:
#### map both RNA-seq reads to the desired genome , treat PE as SE sample 
#	input: 
#		read1 = fastqc.output.trimmedfq1,
#		read2 = fastqc.output.trimmedfq2,
#		readse = fastqc.output.trimmedfqse,
#		index = "{bt2_inx}",
#		thread = "{bt2_thread}"
#	output: 
#		temp("{sample}.step1.bam")
#	run: 
#		if {paired} == "YES" :
#			shell:
#				" {bowtie} --very-sensitive-local -x {input.index} "
#				" --threads {input.thread} "
#				" -U {input.read1}, {input.read2} "
#				" | {samtools} view -bS -F 4 - > {output} "
#		else:
#			shell:
#				" {bowtie} --very-sensitive-local -x {input.index} "
#				" --threads {input.thread} "
#				" -U {input.readse}"
#				" | {samtools} view -bS -F 4 - > {output} "
#


