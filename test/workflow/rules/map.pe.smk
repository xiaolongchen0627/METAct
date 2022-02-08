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