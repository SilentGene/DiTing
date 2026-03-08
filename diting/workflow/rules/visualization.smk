rule visualization:
    input:
        table = os.path.join(out_dir, "pathways_relative_abundance.tab") if not vis_file else vis_file
    output:
        figures = expand("{out_dir}/Figure/{cycle}_{type}.png", 
                   out_dir=[out_dir] if not vis_file else [os.path.dirname(vis_file) or "."], 
                   cycle=["carbon_cycle", "nitrogen_cycle", "DMSP_cycle", "sulfur_cycle"], 
                   type=["sketch"])
    log:
        os.path.join(out_dir, "logs", "visualization.log")
    run:
        from scripts.sketch import sketch
        from scripts.heatmap import heatmap
        import shutil
        
        tab_path = os.path.abspath(input.table)
        tab_dir = os.path.dirname(tab_path)
        tab_name = os.path.basename(tab_path)
        
        cwd = os.getcwd()
        os.chdir(tab_dir)
        
        sketch(tab_name)
        heatmap(tab_name)
        
        fig_dir = "Figure"
        os.makedirs(fig_dir, exist_ok=True)
        os.system(f"mv *.png {fig_dir}/ 2>/dev/null || true")
        os.system(f"mv *.pdf {fig_dir}/ 2>/dev/null || true")
        os.system("rm -rf Figure_tmp heatmap_tmp")
        
        os.chdir(cwd)
