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