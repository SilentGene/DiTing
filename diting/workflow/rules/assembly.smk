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
    threads: threads_cli
    run:
        tmp_dir = os.path.join(out_dir, "assembly_tmp_" + wildcards.sample)
        shell(f"rm -rf {tmp_dir}")
        if use_spades:
            if READS_INTER:
                shell(f"spades.py --meta --12 {input.reads[0]} -t {threads} -o {tmp_dir} -m {memory}")
            else:
                shell(f"spades.py --meta -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -o {tmp_dir} -m {memory}")
            shell(f"cp {tmp_dir}/contigs.fasta {output.fa}")
        else:
            if READS_INTER:
                shell(f"megahit --12 {input.reads[0]} -t {threads} -o {tmp_dir} -f")
            else:
                shell(f"megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -o {tmp_dir} -f")
            shell(f"cp {tmp_dir}/final.contigs.fa {output.fa}")
        shell(f"rm -rf {tmp_dir}")

rule prodigal:
    input:
        assembly = get_assembly
    output:
        faa = os.path.join(out_dir, "ORFs", "{sample}.faa"),
        ffn = os.path.join(out_dir, "ORFs", "{sample}.ffn")
    shell:
        "prodigal -i {input.assembly} -a {output.faa} -d {output.ffn} -p meta -q"
