rule setup_db:
    output:
        ko_list = KO_LIST,
        # We don't strictly "produce" these in the rule anymore as they must exist,
        # but we can use them as triggers or just keep the rule for DMSP parsing.
        ko_profile_sample = os.path.join(KODB_DIR, "K00001.hmm"),
        dmsp_profile = os.path.join(KODB_DIR, "AcuH.hmm")
    log:
        os.path.join(out_dir, "logs", "setup_db.log")
    run:
        from scripts.check import check_kodb, check_DMSP_db
        from scripts.func import DMSP_db_parse
        
        # We assume kofam database exists because the CLI validated it.
        # But we still need to parse DMSP if it's not already in the profiles dir.
        if not check_DMSP_db(KODB_DIR):
            print("Parsing DMSP database...")
            DMSP_db_parse(DMSP_DIR, KODB_DIR, KO_LIST)
