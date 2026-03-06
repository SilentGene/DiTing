rule setup_db:
    output:
        ko_list = os.path.join(KODB_DIR, "ko_list"),
        ko_profile = os.path.join(KODB_DIR, "profiles", "K00001.hmm"),
        dmsp_profile = os.path.join(KODB_DIR, "profiles", "AcuH.hmm")
    run:
        from scripts.check import check_kodb, check_DMSP_db
        from scripts.func import download_db, DMSP_db_parse
        
        if not check_kodb(KODB_DIR):
            print("Downloading KEGG database...")
            download_db(KODB_DIR)
        
        if not check_DMSP_db(KODB_DIR):
            print("Parsing DMSP database...")
            DMSP_db_parse(DMSP_DIR, KODB_DIR)
