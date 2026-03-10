rule kegg_annotation:
    input:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa"),
        ko_list = KO_LIST
    output:
        raw = temp(os.path.join(out_dir, "KEGG_annotation", "pfamscan_out", "{sample}-ko-annotations.tsv")),
        filtered = os.path.join(out_dir, "KEGG_annotation", "pfamscan_out", "{sample}-ko-annotations-filtered.tsv")
    threads: max(2, min(8, int(threads_cli) // max(1, len(BASENAMES))))
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
        filtered = os.path.join(out_dir, "DMSP_annotation", "hmmout", "{sample}-dmsp-annotations-filtered.tsv")
    threads: max(2, min(8, int(threads_cli) // max(1, len(BASENAMES))))
    log:
        os.path.join(out_dir, "logs", "dmsp_annotation_{sample}.log")
    run:
        import os
        hmmout_dir = os.path.join(out_dir, "DMSP_annotation", "hmmout")
        os.makedirs(hmmout_dir, exist_ok=True)
        dmsp_profiles_dir = os.path.join(DMSP_DIR, "profiles")
        dmsp_list = os.path.join(DMSP_DIR, "DMSP_related_gene.list")
        
        # Read evalues from dmsp_list
        evalues = {}
        with open(dmsp_list, 'r') as f:
            for line in f:
                if line.startswith('knum'): continue
                if not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    knum = parts[0]
                    evalue = parts[1]
                    evalues[knum] = evalue
        
        # Run hmmsearch per profile
        for knum, evalue in evalues.items():
            hmm_profile = os.path.join(dmsp_profiles_dir, f"{knum}.hmm")
            if os.path.exists(hmm_profile):
                out_hmmout = os.path.join(hmmout_dir, f"{knum}.{wildcards.sample}.hmmout")
                shell(f"hmmsearch -o /dev/null --tblout {out_hmmout} -E {evalue} --cpu {threads} {hmm_profile} {input.faa} >> {{log}} 2>&1")
        
        # Filter best hits
        filter_script = os.path.join(workflow.basedir, "scripts", "dmsp_filter.py")
        shell(f"python {filter_script} -i {hmmout_dir} -o {output.filtered} -s {wildcards.sample} >> {{log}} 2>&1")
