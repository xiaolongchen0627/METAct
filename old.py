



#else :
#    rule trim_galore:
#        input:
#            trim_inputs
#        output:
#            temp(os.path.join(config["results_dir"],"{sample}",'trimmed',"{sample}_trimmed.fq.gz"))
#        params:
#            basename = "{sample}",
#            trimmed_dir = os.path.join(config["results_dir"],"{sample}","trimmed")
#        shell: 
#            "trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} {input} ; " 
#            "touch test.finish"       
##else:
#   if is_paired:
#       rule ln_fq:
#           input:
#               trim_inputs
#           output:
#               temp(os.path.join(config["results_dir"],"{sample}",'trimmed',"{sample}_val_1.fq.gz"))
#               temp(os.path.join(config["results_dir"],"{sample}",'trimmed',"{sample}_val_2.fq.gz"))
#           params:
#               basename = "{sample}",
#               outputDir = os.path.join(config["results_dir"],"{sample}","trimmed")
#           shell:
#                


#rule all : 
#    input:
 #       "test.finish"
        #expand("{sample}/{sample}.step1.sam",sample=samples["sample_name"])
       

#def is_paired(sample):
#    sp = samples.loc[sample]
#    return sp["fq1"]
#

#def need_trim(sample):
#    sp = samples.loc[sample]
#    return sp["need_trim"]

#def get_fq(wildcards):
#    if need_trim(wildcards.sample):
#        if is_paired(wildcards.sample):
#            output.fq = dict(
#                zip(
#                    ["fq1","fq2"],
#                    expand("{wd}/{sample}/{sample}_val_{group}.fq.gz",
#                        group=['1','2'],
#                        **wildcards)
#                )
#            )
#            input.fq = {"fq1":wildcards.sample["fq1"],
#                    "fq2":wildcards.sample["fq2"]}
#            
#                    
#        else:
#            output.fq = dict(
#                zip(
#                    ["fq1"],
#                    expand("{wd}/{sample}/{sample}_trimmed.fq.gz".format(**wildcards))
#                )
#            )
#            input.fq ={"fq1":wildcards.sample["fq1"]}
#        print(input.fq)
#        print(output.fq)
#    return input.fq


#rule all : 
#    input:
#        expand("{wd}/{sample}/{sample}.step1.sam",wd=config["wd"],sample=samples["sample_name"])
#       
#rule test:
#    input:
#        ""
#    output:
#        "test.txt"
#    shell:
#        "touch test.txt"
#

#rule test:
#    input:
#        unpack(get_fq)
#    output:
#        "test.txt"
#    shell:
#        "touch test.txt"
#if paired :
#    rule TrimGalore:
#        input:
#            r1 = {read1},
#            r2 = {read2}
#        output:
#            r1 = temp("{wd}/{sample}/{sample}_val_1.fq.gz"),
#            r2 = temp("{wd}/{sample}/{sample}_val_2.fq.gz")
#        params:
#            basename = "{sample}",
#            outputDir = "{wd}/{sample}"
#        shell: 
#            "trim_galore --gzip --output_dir {params.outputDir} --basename {params.basename} --paired {input.r1} {input.r2} "
#          
#    rule Alignment:
#        input:
#            r1 = rules.TrimGalore.output.r1,
#            r2 = rules.TrimGalore.output.r2
#        output:
#            bam = temp("{wd}/{sample}/{sample}.step1.sam")
#        shell : 
#            " bowtie2 --very-sensitive-local"
#            " --thread {bt2_thread}"
#            " -x {bt2_idx} -U {input.r1},{input.r2} "
#            " > {output} "
#
#else :
#    rule TrimGalore:
#        input:
#            r1 = "{read1}"
#        output:
#            r1 = temp("{wd}/{sample}/{sample}_trimmed.fq.gz")
#        params:
#            basename = "{sample}",
#            outputDir = "{wd}/{sample}"
#        shell: 
#            "trim_galore --output_dir {params.outputDir} --basename {params.basename} {input.r1}"
#
#    rule Alignment:
#        input:
#            r1 = rules.TrimGalore.output.r1
#        output:
#            bam = temp("{wd}/{sample}/{sample}.step1.sam")
#        shell : 
#            " bowtie2 --very-sensitive-local"
#            " --thread {bt2_thread}"
#            " -x {bt2_idx} -U {input.r1}"
#            " > {output} "
#
#rule filterSam:
#    input:
#        rules.Alignment.output
#    output:
#        bam = temp("{wd}/{sample}/{sample}.step1.filter.bam")
#rule bam2Fastq:
#    input:
#        rules.Alignment.output
#    output:
#        fqStep2 = temp("{wd}/{sample}/{sample}.step2.fq")
#    shell:
#        "picard SamToFastq "
#        "I={input} "
#        "F={output} "
#        "VALIDATION_STRINGENCY=LENIENT "

