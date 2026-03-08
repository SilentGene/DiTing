rule merge_annotations:
    input:
        kegg_filtered = expand(os.path.join(out_dir, "KEGG_annotation", "pfamscan_out", "{sample}-ko-annotations-filtered.tsv"), sample=BASENAMES),
        dmsp_filtered = expand(os.path.join(out_dir, "DMSP_annotation", "hmmout", "{sample}-dmsp-annotations-filtered.tsv"), sample=BASENAMES)
    output:
        ko_merged = os.path.join(out_dir, "KEGG_annotation", "ko_merged.txt")
    log:
        os.path.join(out_dir, "logs", "merge_annotations.log")
    run:
        with open(output.ko_merged, 'w') as fo:
            fo.write('#sample\tgene_id\tk_number\n')
            for kegg_f, dmsp_f in zip(input.kegg_filtered, input.dmsp_filtered):
                kegg_hits = {}
                with open(kegg_f, 'r') as fi:
                    for line in fi:
                        if line.startswith('#'): continue
                        if not line.strip(): continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            gene = parts[1]
                            kegg_hits[gene] = line
                
                dmsp_hits = {}
                with open(dmsp_f, 'r') as fi:
                    for line in fi:
                        if line.startswith('#'): continue
                        if not line.strip(): continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            gene = parts[1]
                            dmsp_hits[gene] = line
                
                all_genes = sorted(list(set(kegg_hits.keys()) | set(dmsp_hits.keys())))
                for gene in all_genes:
                    if gene in dmsp_hits:
                        fo.write(dmsp_hits[gene])
                    else:
                        fo.write(kegg_hits[gene])

rule merge_abun_ko:
    input:
        abuns = expand(os.path.join(out_dir, "Abundance", "{sample}.abundance"), sample=BASENAMES),
        ko_merged = rules.merge_annotations.output.ko_merged
    output:
        ko_abun = os.path.join(out_dir, "KEGG_annotation", "ko_abun.txt")
    log:
        os.path.join(out_dir, "logs", "merge_abun_ko.log")
    run:
        import os
        from scripts.func import merge_abun_ko
        merge_abun_ko(os.path.join(out_dir, "Abundance"), input.ko_merged, output.ko_abun)

rule process_tables:
    input:
        ko_abun = rules.merge_abun_ko.output.ko_abun,
        faas = expand(os.path.join(out_dir, "ORFs", "{sample}.faa"), sample=BASENAMES)
    output:
        ko_abundance_among_samples = os.path.join(out_dir, "ko_abundance_among_samples.tab"),
        pathways_table = os.path.join(out_dir, "pathways_relative_abundance.tab")
    log:
        os.path.join(out_dir, "logs", "process_tables.log")
    run:
        import os
        from scripts.func import table_of_ko_abundance_among_samples, build_gene_family, hierarchical_ko_abundance_among_samples, pathway_parser
        from scripts.transposition import transposition
        
        ko_abundance_among_samples = os.path.join(out_dir, "KEGG_annotation", "ko_abundance_among_samples.tab")
        table_of_ko_abundance_among_samples(input.ko_abun, ko_abundance_among_samples)
        
        family_dir = os.path.join(out_dir, "Gene_family")
        os.makedirs(family_dir, exist_ok=True)
        build_gene_family(os.path.join(out_dir, "ORFs"), input.ko_abun, family_dir)
        
        hier_tab = os.path.join(out_dir, "pathways_relative_abundance_gene_level.tab")
        cycle_tab = os.path.join(TABLE_DIR, "KO_affilated_to_biogeochemical_cycle.tab")
        hierarchical_ko_abundance_among_samples(ko_abundance_among_samples, cycle_tab, hier_tab)
        
        pathway_tmp = os.path.join(out_dir, "pathways_relative_abundance_tmp.tab")
        pathway_parser(input.ko_abun, pathway_tmp)
        transposition(pathway_tmp)
        os.rename(pathway_tmp + ".transposition", output.pathways_table)
        os.rename(ko_abundance_among_samples, output.ko_abundance_among_samples)
        os.remove(pathway_tmp)
