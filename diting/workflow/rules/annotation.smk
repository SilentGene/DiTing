rule kegg_annotation:
    input:
        faa = NR_FAA,
        ko_list = lambda wildcards: KO_LIST
    output:
        raw = temp(os.path.join(
            out_dir, "KEGG_annotation", "pfamscan_out",
            "catalogue-ko-annotations.tsv"
        )),
        filtered = os.path.join(
            out_dir, "KEGG_annotation", "pfamscan_out",
            "catalogue-ko-annotations-filtered.tsv"
        )
    threads: threads_cli
    log:
        os.path.join(out_dir, "logs", "kegg_annotation_catalogue.log")
    run:
        import os
        import shlex

        pfamscan_out_dir = os.path.join(out_dir, "KEGG_annotation", "pfamscan_out")
        os.makedirs(pfamscan_out_dir, exist_ok=True)
        tmp_dir = os.path.join(pfamscan_out_dir, "catalogue-ko-tmp")
        command = (
            f"exec_annotation -p {shlex.quote(str(KODB_DIR))} "
            f"-k {shlex.quote(str(input.ko_list))} --cpu {threads} "
            "-f detail-tsv --e-value 1e-5 "
            f"--tmp-dir {shlex.quote(tmp_dir)} "
            f"-o {shlex.quote(str(output.raw))} {shlex.quote(str(input.faa))} "
            f"> {shlex.quote(str(log))} 2>&1"
        )
        shell(command)

        filter_script = os.path.join(workflow.basedir, "scripts", "kofamscan_filter.py")
        command = (
            f"python {shlex.quote(filter_script)} "
            f"-i {shlex.quote(str(output.raw))} "
            f"-o {shlex.quote(str(output.filtered))} "
            f"-E 1e-5 >> {shlex.quote(str(log))} 2>&1"
        )
        shell(command)
        shell(f"rm -rf {shlex.quote(tmp_dir)}")


rule dmsp_annotation:
    input:
        faa = NR_FAA
    output:
        filtered = os.path.join(
            out_dir, "DMSP_annotation", "hmmout",
            "catalogue-dmsp-annotations-filtered.tsv"
        )
    threads: threads_cli
    log:
        os.path.join(out_dir, "logs", "dmsp_annotation_catalogue.log")
    run:
        import os
        import shlex

        hmmout_dir = os.path.join(out_dir, "DMSP_annotation", "hmmout", "catalogue")
        os.makedirs(hmmout_dir, exist_ok=True)
        dmsp_profiles_dir = os.path.join(DMSP_DIR, "profiles")
        dmsp_list = os.path.join(DMSP_DIR, "DMSP_related_gene.list")

        evalues = {}
        with open(dmsp_list, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("knum") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    evalues[parts[0]] = parts[1]

        for knum, evalue in evalues.items():
            hmm_profile = os.path.join(dmsp_profiles_dir, f"{knum}.hmm")
            if os.path.exists(hmm_profile):
                hmmout = os.path.join(hmmout_dir, f"{knum}.catalogue.hmmout")
                command = (
                    "hmmsearch -o /dev/null "
                    f"--tblout {shlex.quote(hmmout)} -E {evalue} "
                    f"--cpu {threads} {shlex.quote(hmm_profile)} "
                    f"{shlex.quote(str(input.faa))} "
                    f">> {shlex.quote(str(log))} 2>&1"
                )
                shell(command)

        filter_script = os.path.join(workflow.basedir, "scripts", "dmsp_filter.py")
        command = (
            f"python {shlex.quote(filter_script)} -i {shlex.quote(hmmout_dir)} "
            f"-o {shlex.quote(str(output.filtered))} "
            f">> {shlex.quote(str(log))} 2>&1"
        )
        shell(command)
