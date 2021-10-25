rule all : 
    input:
        "{wd}/{sample}/{sample}.step1.bam":

if paired == "PE" :
    rule TrimGalore:
        input:
            r1 = "{read1}",
            r2 = "{read2}"
        output:
            r1 = temp("{wd}/{sample}/{sample}_val_1.fq.gz"),
            r2 = temp("{wd}/{sample}/{sample}_val_2.fq.gz")
        params:
            basename = "{sample}",
            outputDir = "{wd}/{sample}"
        log:
            out = "{wd}/logs/TrimGalore.{sample}.out",
            err = "{wd}/logs/TrimGalore.{sample}.err"
        shell: 
            "trim_galore --gzip --output_dir {params.outputDir} --basename {params.basename} --paired {input.r1} {input.r2} > {log.out} 2> {log.err}"
          
    rule Alignment:
        input:
            r1 = "{TrimGalore.output.r1}",
            r2 = "{TrimGalore.output.r2}"
        output:
            bam = "{wd}/{sample}/{sample}.step1.bam"
        shell : 
            " bowtie2 --very-sensitive-local"
            " --thread {bt2_thread}"
            " -x {bt2_idx} -U {input.r1}, {input.r2} "
            " | samtools view -bS -F 4 - >{output} "

else :
    rule TrimGalore:
        input:
            r1 = "{read1}"
        output:
            r1 = temp("{wd}/{sample}/{sample}_trimmed.fq.gz")
        params:
            basename = "{sample}",
            outputDir = "{wd}/{sample}"
        log:
            out = "{wd}/logs/TrimGalore.{sample}.out",
            err = "{wd}/logs/TrimGalore.{sample}.err"
        shell: 
            "trim_galore --output_dir {params.outputDir} --basename {params.basename} {input.r1}  > {log.out} 2> {log.err}"

    rule Alignment:
        input:
            r1 = "{TrimGalore.output.r1}"
        output:
            bam = "{wd}/{sample}/{sample}.step1.bam"
        shell : 
            " bowtie2 --very-sensitive-local"
            " --thread {bt2_thread}"
            " -x {bt2_idx} -U {input.r1}"
            " | samtools view -bS -F 4 - >{output} "