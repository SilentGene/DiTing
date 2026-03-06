rule bwa_index:
    input:
        ffn = os.path.join(out_dir, "ORFs", "{sample}.ffn")
    output:
        multiext(os.path.join(out_dir, "BBMap", "bwa_index", "{sample}"), ".amb", ".ann", ".bwt", ".pac", ".sa")
    params:
        prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
    shell:
        "bwa index -p {params.prefix} {input.ffn} > /dev/null 2>&1"

rule bwa_mem:
    input:
        reads = get_reads,
        idx = rules.bwa_index.output,
        ffn = os.path.join(out_dir, "ORFs", "{sample}.ffn")
    output:
        sam = temp(os.path.join(out_dir, "BBMap", "mapping", "{sample}.sam")) if not noclean else os.path.join(out_dir, "BBMap", "mapping", "{sample}.sam")
    threads: threads_cli
    params:
        prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0]
    run:
        if READS_INTER:
            shell(f"bwa mem -p -t {threads} {{params.prefix}} {{input.reads[0]}} > {{output.sam}}")
        else:
            shell(f"bwa mem -t {threads} {{params.prefix}} {{input.reads[0]}} {{input.reads[1]}} > {{output.sam}}")

rule bbmap_pileup:
    input:
        sam = rules.bwa_mem.output.sam
    output:
        pileup = os.path.join(out_dir, "BBMap", "coverage", "{sample}.pileup")
    shell:
        "pileup.sh in={input.sam} out={output.pileup} -Xmx20g > /dev/null 2>&1"

rule calculate_abundance:
    input:
        pileup = rules.bbmap_pileup.output.pileup
    output:
        abun = os.path.join(out_dir, "Abundance", "{sample}.abundance")
    run:
        from scripts.func import gene_relative_abun
        gene_relative_abun(input.pileup, wildcards.sample, os.path.dirname(output.abun))
