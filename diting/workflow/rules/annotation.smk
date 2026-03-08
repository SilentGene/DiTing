rule kegg_annotation:
    input:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa"),
        ko_list = KO_LIST
    output:
        raw = temp(os.path.join(out_dir, "KEGG_annotation", "pfamscan_out", "{sample}-ko-annotations.tsv")),
        filtered = os.path.join(out_dir, "KEGG_annotation", "pfamscan_out", "{sample}-ko-annotations-filtered.tsv")
    threads: 4 
    log:
        os.path.join(out_dir, "logs", "kegg_annotation_{sample}.log")
    run:
        import os
        pfamscan_out_dir = os.path.join(out_dir, "KEGG_annotation", "pfamscan_out")
        os.makedirs(pfamscan_out_dir, exist_ok=True)
        tmp_dir = os.path.join(pfamscan_out_dir, f"{wildcards.sample}-ko-tmp")
        profiles_dir = KODB_DIR
        
        shell(f"exec_annotation -p {profiles_dir} -k {input.ko_list} --cpu {threads} -f detail-tsv --e-value 1e-5 --tmp-dir {tmp_dir} -o {output.raw} {input.faa} > {{log}} 2>&1")
        
        filter_script = os.path.join(workflow.basedir, "scripts", "kofamscan_filter.py")
        shell(f"python {filter_script} -i {output.raw} -o {output.filtered} -s {wildcards.sample} -E 1e-5 >> {{log}} 2>&1")
        
        shell(f"rm -rf {tmp_dir}")

rule dmsp_annotation:
    input:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa")
    output:
        raw = temp(os.path.join(out_dir, "DMSP_annotation", "pfamscan_out", "{sample}-dmsp-annotations.tsv")),
        filtered = os.path.join(out_dir, "DMSP_annotation", "pfamscan_out", "{sample}-dmsp-annotations-filtered.tsv")
    threads: 4 
    log:
        os.path.join(out_dir, "logs", "dmsp_annotation_{sample}.log")
    run:
        import os
        pfamscan_out_dir = os.path.join(out_dir, "DMSP_annotation", "pfamscan_out")
        os.makedirs(pfamscan_out_dir, exist_ok=True)
        tmp_dir = os.path.join(pfamscan_out_dir, f"{wildcards.sample}-dmsp-tmp")
        dmsp_profiles_dir = os.path.join(DMSP_DIR, "profiles")
        dmsp_list = os.path.join(DMSP_DIR, "DMSP_related_gene.list")
        
        shell(f"exec_annotation -p {dmsp_profiles_dir} -k {dmsp_list} --cpu {threads} -f detail-tsv --e-value 1e-5 --tmp-dir {tmp_dir} -o {output.raw} {input.faa} > {{log}} 2>&1")
        
        filter_script = os.path.join(workflow.basedir, "scripts", "kofamscan_filter.py")
        shell(f"python {filter_script} -i {output.raw} -o {output.filtered} -s {wildcards.sample} -E 1e-5 >> {{log}} 2>&1")
        
        shell(f"rm -rf {tmp_dir}")
