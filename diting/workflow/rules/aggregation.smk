rule merge_annotations:
    input:
        kegg_filtered = rules.kegg_annotation.output.filtered,
        dmsp_filtered = rules.dmsp_annotation.output.filtered
    output:
        ko_merged = os.path.join(out_dir, "KEGG_annotation", "ko_merged.txt")
    log:
        os.path.join(out_dir, "logs", "merge_annotations.log")
    run:
        def read_annotations(path):
            annotations = {}
            with open(path, encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        annotations[parts[0]] = parts[1]
            return annotations

        kegg_hits = read_annotations(input.kegg_filtered)
        dmsp_hits = read_annotations(input.dmsp_filtered)
        all_genes = sorted(set(kegg_hits) | set(dmsp_hits))

        with open(output.ko_merged, "w", encoding="utf-8") as output_handle:
            output_handle.write("#gene_id\tk_number\n")
            for gene_id in all_genes:
                k_number = dmsp_hits.get(gene_id, kegg_hits.get(gene_id))
                output_handle.write(f"{gene_id}\t{k_number}\n")


rule merge_abun_ko:
    input:
        matrix = rules.normalize_coverm_tpm.output.matrix,
        ko_merged = rules.merge_annotations.output.ko_merged
    output:
        ko_abun = os.path.join(out_dir, "KEGG_annotation", "ko_abun.txt")
    log:
        os.path.join(out_dir, "logs", "merge_abun_ko.log")
    run:
        from scripts.func import merge_abun_ko

        merge_abun_ko(input.matrix, input.ko_merged, output.ko_abun)


rule process_tables:
    input:
        ko_abun = rules.merge_abun_ko.output.ko_abun,
        ko_merged = rules.merge_annotations.output.ko_merged,
        faa = NR_FAA
    output:
        ko_abundance_among_samples = os.path.join(
            out_dir, "ko_abundance_among_samples.tab"
        ),
        pathways_table = os.path.join(out_dir, "pathways_relative_abundance.tab")
    log:
        os.path.join(out_dir, "logs", "process_tables.log")
    run:
        import os

        from scripts.func import (
            build_gene_family,
            hierarchical_ko_abundance_among_samples,
            pathway_parser,
            table_of_ko_abundance_among_samples,
        )
        from scripts.transposition import transposition

        intermediate_ko_table = os.path.join(
            out_dir, "KEGG_annotation", "ko_abundance_among_samples.tab"
        )
        table_of_ko_abundance_among_samples(input.ko_abun, intermediate_ko_table)

        family_dir = os.path.join(out_dir, "Gene_family")
        os.makedirs(family_dir, exist_ok=True)
        build_gene_family(input.faa, input.ko_merged, family_dir)

        hierarchy_table = os.path.join(
            out_dir, "pathways_relative_abundance_gene_level.tab"
        )
        cycle_table = os.path.join(
            TABLE_DIR, "KO_affilated_to_biogeochemical_cycle.tab"
        )
        hierarchical_ko_abundance_among_samples(
            intermediate_ko_table, cycle_table, hierarchy_table
        )

        pathway_tmp = os.path.join(out_dir, "pathways_relative_abundance_tmp.tab")
        pathway_parser(input.ko_abun, pathway_tmp)
        transposition(pathway_tmp)
        os.replace(pathway_tmp + ".transposition", output.pathways_table)
        os.replace(intermediate_ko_table, output.ko_abundance_among_samples)
        os.remove(pathway_tmp)
