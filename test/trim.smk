if is_paired:
    rule trim_galore:
        input:
            trim_inputs
        output:
            os.path.join(config["working_dir"],"trimmed","{sample}_val_1.fq.gz"),
            os.path.join(config["working_dir"],"trimmed","{sample}_val_2.fq.gz")
        params:
            basename = "{sample}",
            trimmed_dir = config["working_dir"]+"/trimmed/",
            tmpdir = config["working_dir"]+"/trimmed/"
        shell: 
            "trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} --paired {input} ; "
            "echo trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} --paired {input} ; "
    
    rule align_step1:
        input:
            fq1 = config["working_dir"]+"/trimmed/"+"{sample}"+"_val_1.fq.gz",
            fq2 = config["working_dir"]+"/trimmed/"+"{sample}"+"_val_2.fq.gz"
        output:
            os.path.join(config["working_dir"],"trimmed","{sample}.sam")
        params:
            bt2_thread = config['bt2_thread'],
            bt2_idx = config['bt2_idx']
        shell:
            " bowtie2 --very-sensitive-local"
            " --thread {params.bt2_thread}"
            " -x {params.bt2_idx} -U {input.fq1},{input.fq2} -S "
            " {output} "
else:
    rule trim_galore:
        input:
            trim_inputs
        output:
            os.path.join(config["working_dir"],"trimmed","{sample}_trimmed.fq.gz")
        params:
            basename = "{sample}",
            trimmed_dir = config["working_dir"]+"/trimmed/",
            tmpdir = config["working_dir"]+"/trimmed/"
        shell: 
            "trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} {input} ; "
            "echo trim_galore --gzip --output_dir {params.trimmed_dir} --basename {params.basename} {input} ; "
    rule align_step1:
        input:
            trim_outputs
        output:
            os.path.join(config["working_dir"],"trimmed","{sample}.sam")
        params:
            bt2_thread = config['bt2_thread'],
            bt2_idx = config['bt2_idx']
        shell:
            " bowtie2 --very-sensitive-local"
            " --thread {params.bt2_thread}"
            " -x {params.bt2_idx} -U {input} -S "
            " {output} "  
