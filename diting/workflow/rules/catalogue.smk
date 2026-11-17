rule merge_sample_orfs:
    input:
        ffns = expand(os.path.join(out_dir, "ORFs", "{sample}.ffn"), sample=BASENAMES),
        faas = expand(os.path.join(out_dir, "ORFs", "{sample}.faa"), sample=BASENAMES)
    output:
        ffn = os.path.join(NR_DIR, "all_genes.ffn"),
        faa = os.path.join(NR_DIR, "all_genes.faa"),
        id_map = os.path.join(NR_DIR, "gene_id_map.tsv")
    log:
        os.path.join(out_dir, "logs", "merge_sample_orfs.log")
    params:
        script = os.path.join(workflow.basedir, "scripts", "merge_orfs.py"),
        orf_dir = os.path.join(out_dir, "ORFs"),
        samples = BASENAMES
    shell:
        """
        python {params.script:q} --orf-dir {params.orf_dir:q} --samples {params.samples:q} --output-ffn {output.ffn:q} --output-faa {output.faa:q} --id-map {output.id_map:q} > {log:q} 2>&1
        """


rule filter_short_orfs:
    input:
        ffn = rules.merge_sample_orfs.output.ffn
    output:
        ffn = os.path.join(NR_DIR, "filtered_100bp.ffn")
    log:
        os.path.join(out_dir, "logs", "filter_short_orfs.log")
    params:
        script = os.path.join(workflow.basedir, "scripts", "filter_short_orfs.py")
    shell:
        """
        python {params.script:q} --input {input.ffn:q} --output {output.ffn:q} --minimum-length 100 > {log:q} 2>&1
        """


rule mmseqs_linclust:
    input:
        ffn = rules.filter_short_orfs.output.ffn
    output:
        ffn = NR_FFN,
        clusters = os.path.join(NR_DIR, f"{NR_PREFIX}_cluster.tsv")
    threads: threads_cli
    log:
        os.path.join(out_dir, "logs", "mmseqs_linclust.log")
    params:
        prefix = os.path.join(NR_DIR, NR_PREFIX),
        tmp = os.path.join(NR_DIR, f"{NR_PREFIX}.tmp"),
        identity = min_seq_id
    shell:
        """
        mmseqs easy-linclust {input.ffn:q} {params.prefix:q} {params.tmp:q} --min-seq-id {params.identity} --cov-mode 1 -c 0.9 -a --threads {threads} > {log:q} 2>&1
        mv {params.prefix}_rep_seq.fasta {output.ffn:q}
        rm -rf {params.tmp:q}
        """


rule catalogue_representative_names:
    input:
        ffn = rules.mmseqs_linclust.output.ffn
    output:
        names = os.path.join(NR_DIR, f"{NR_PREFIX}_100bp.ffn.names")
    run:
        from scripts.fasta_io import read_fasta

        with open(output.names, "w", encoding="utf-8") as handle:
            for record in read_fasta(input.ffn):
                handle.write(record.identifier + "\n")


rule index_all_catalogue_proteins:
    input:
        faa = rules.merge_sample_orfs.output.faa
    output:
        ssi = os.path.join(NR_DIR, "all_genes.faa.ssi")
    log:
        os.path.join(out_dir, "logs", "index_all_catalogue_proteins.log")
    shell:
        """
        rm -f {output.ssi:q}
        esl-sfetch --index {input.faa:q} > {log:q} 2>&1
        """


rule extract_catalogue_proteins:
    input:
        faa = rules.merge_sample_orfs.output.faa,
        ssi = rules.index_all_catalogue_proteins.output.ssi,
        names = rules.catalogue_representative_names.output.names
    output:
        faa = NR_FAA
    log:
        os.path.join(out_dir, "logs", "extract_catalogue_proteins.log")
    shell:
        "esl-sfetch -f {input.faa:q} {input.names:q} > {output.faa:q} 2> {log:q}"
