def get_assembly(wildcards):
    if assembly_dir:
        return os.path.join(assembly_dir, wildcards.sample + ".fa")
    return os.path.join(out_dir, "Assembly", wildcards.sample + ".fa")

def get_reads(wildcards):
    if READS_INTER:
        return [os.path.join(reads_dir, wildcards.sample + READS_SUF)]
    return [os.path.join(reads_dir, wildcards.sample + "_1" + READS_SUF),
            os.path.join(reads_dir, wildcards.sample + "_2" + READS_SUF)]

rule assemble:
    input:
        reads = get_reads
    output:
        fa = os.path.join(out_dir, "Assembly", "{sample}.fa")
    threads: 8
    log:
        os.path.join(out_dir, "logs", "assemble_{sample}.log")
    run:
        tmp_dir = os.path.join(out_dir, "assembly_tmp_" + wildcards.sample)
        shell(f"rm -rf {tmp_dir}")
        if use_spades:
            if READS_INTER:
                shell(f"spades.py --meta --12 {input.reads[0]} -t {threads} -o {tmp_dir} -m {memory} > {{log}} 2>&1")
            else:
                shell(f"spades.py --meta -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -o {tmp_dir} -m {memory} > {{log}} 2>&1")
            shell(f"cp {tmp_dir}/contigs.fasta {output.fa} >> {{log}} 2>&1")
        else:
            if READS_INTER:
                shell(f"megahit --12 {input.reads[0]} -t {threads} -o {tmp_dir} -f > {{log}} 2>&1")
            else:
                shell(f"megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -o {tmp_dir} -f > {{log}} 2>&1")
            shell(f"cp {tmp_dir}/final.contigs.fa {output.fa} >> {{log}} 2>&1")
        shell(f"rm -rf {tmp_dir}")

rule prodigal:
    input:
        assembly = get_assembly
    output:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa"),
        ffn = os.path.join(out_dir, "ORFs", "{sample}.ffn")
    log:
        os.path.join(out_dir, "logs", "prodigal_{sample}.log")
    shell:
        "prodigal -i {input.assembly} -a {output.faa} -d {output.ffn} -p meta -q > {log} 2>&1"
