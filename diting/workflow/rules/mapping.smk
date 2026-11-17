def get_all_reads():
    reads = []
    for sample in BASENAMES:
        reads.extend(get_reads(type("Wildcards", (), {"sample": sample})()))
    return reads


rule coverm_tpm:
    input:
        ffn = NR_FFN,
        reads = get_all_reads()
    output:
        raw = temp(os.path.join(out_dir, "Abundance", "coverm_tpm.raw.tsv"))
    threads: threads_cli
    log:
        os.path.join(out_dir, "logs", "coverm_tpm.log")
    params:
        reads = get_all_reads(),
        bam_cache = os.path.join(out_dir, "BBMap", "mapping")
    run:
        import os
        import shlex

        read_option = "--interleaved" if READS_INTER else "--coupled"
        quoted_reads = " ".join(shlex.quote(str(path)) for path in params.reads)
        cache_option = ""
        if noclean:
            os.makedirs(params.bam_cache, exist_ok=True)
            cache_option = (
                "--bam-file-cache-directory "
                + shlex.quote(str(params.bam_cache))
            )
        command = (
            f"coverm contig {read_option} {quoted_reads} "
            f"--reference {shlex.quote(str(input.ffn))} --threads {threads} "
            "--methods tpm --min-mapq 1 --output-format dense "
            f"--output-file {shlex.quote(str(output.raw))} {cache_option} "
            f"> {shlex.quote(str(log))} 2>&1"
        )
        shell(command)


rule normalize_coverm_tpm:
    input:
        raw = rules.coverm_tpm.output.raw,
        ffn = NR_FFN
    output:
        matrix = os.path.join(out_dir, "gene_abundance_among_samples.tab"),
        abundances = expand(
            os.path.join(out_dir, "Abundance", "{sample}.abundance"),
            sample=BASENAMES
        )
    log:
        os.path.join(out_dir, "logs", "normalize_coverm_tpm.log")
    params:
        script = os.path.join(workflow.basedir, "scripts", "normalize_coverm_tpm.py"),
        samples = BASENAMES,
        abundance_dir = os.path.join(out_dir, "Abundance")
    shell:
        """
        python {params.script:q} --input {input.raw:q} --catalogue {input.ffn:q} --samples {params.samples:q} --output-matrix {output.matrix:q} --abundance-dir {params.abundance_dir:q} > {log:q} 2>&1
        """
