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