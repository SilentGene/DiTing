rule kegg_annotation:
    input:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa"),
        ko_list = KO_LIST
    output:
        raw = temp(os.path.join(out_dir, "KEGG_annotation", "hmmout", "{sample}-ko-annotations.tsv")),
        filtered = os.path.join(out_dir, "KEGG_annotation", "hmmout", "{sample}-ko-annotations-filtered.tsv")
    threads: threads_cli   
    run:
        hmmout_dir = os.path.join(out_dir, "KEGG_annotation", "hmmout")
        os.makedirs(hmmout_dir, exist_ok=True)
        tmp_dir = os.path.join(hmmout_dir, f"{wildcards.sample}-ko-tmp")
        profiles_dir = KODB_DIR
        
        shell(f"exec_annotation -p {profiles_dir} -k {input.ko_list} --cpu {threads} -f detail-tsv --e-value 1e-5 --tmp-dir {tmp_dir} -o {output.raw} {input.faa}")
        
        filter_script = os.path.join(workflow.basedir, "scripts", "kofamscan_filter.py")
        shell(f"python {filter_script} -i {output.raw} -o {output.filtered} -s {wildcards.sample} -E 1e-5")
        
        shell(f"rm -rf {tmp_dir}")
