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