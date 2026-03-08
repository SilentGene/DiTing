"""
Define variables
"""

import os
import sys
from .args import *

# Defensive access to args
arg_o = getattr(args, 'o', None)
arg_r = getattr(args, 'r', None)
arg_n = getattr(args, 'n', 4)
arg_m = getattr(args, 'm', 50)
arg_a = getattr(args, 'a', None)
arg_vis = getattr(args, 'vis', False)

if arg_vis:
    ABUNDANCE_TABLE = arg_vis
else:
    READS_DIR = arg_r
    OUT_DIR = arg_o
    THREADS = arg_n
    MEM = arg_m
    
    if OUT_DIR:
        ASSEMBLY_TMP = os.path.join(OUT_DIR, 'assembly_tmp')
        PRODIGAL_DIR = os.path.join(OUT_DIR, 'ORFs')
        BBMAP_DIR = os.path.join(OUT_DIR, 'BBMap')
        GENE_ABUN_DIR = os.path.join(OUT_DIR, 'Abundance')
        KEGG_DIR = os.path.join(OUT_DIR, 'KEGG_annotation')
        GENE_FAMILY = os.path.join(OUT_DIR, 'Gene_family')
        ASSEMBLY_DIR = arg_a if arg_a else os.path.join(OUT_DIR, 'Assembly')
    else:
        ASSEMBLY_TMP = None
        PRODIGAL_DIR = None
        BBMAP_DIR = None
        GENE_ABUN_DIR = None
        KEGG_DIR = None
        GENE_FAMILY = None
        ASSEMBLY_DIR = arg_a

    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    KODB_DIR = os.path.join(ROOT_DIR, 'kofam_database') if ROOT_DIR else 'kofam_database'
    DMSP_DIR = os.path.join(ROOT_DIR, 'DMSP_database') if ROOT_DIR else 'DMSP_database'
    TABLE = os.path.join(ROOT_DIR, 'table') if ROOT_DIR else 'table'
    READS_INTER = None

    BASENAMES = []
    READS_SUF = ''
    ASSEMBLY_SUF = 'fa'
    ABUNDANCE_TABLE = 'pathways_relative_abundance.tab'