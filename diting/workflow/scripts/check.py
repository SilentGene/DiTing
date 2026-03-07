"""
Check data integrity
"""
import os
import re
from scripts.logformatter import *


def check_kodb(ko_db):
    """
    Check if the kofam database is valid.
    ko_db is now the directory containing .hmm files.
    """
    if not os.path.isdir(ko_db):
        return False
    # Check if there is at least one .hmm file
    import glob
    if not glob.glob(os.path.join(ko_db, "*.hmm")):
        return False
    return True

def check_DMSP_db(ko_db):
    hmm_path = os.path.join(ko_db, 'AcuH.hmm')
    if not os.path.exists(hmm_path):
        return False
    return True


def check_reads_assembly(assembly_dir, basenames):
    """
    Check if the assembly provided can match corresponding reads pairs

    :param str assembly_dir: directory containing metagenomic assemblies in fasta format
    :param list basenames: a list of basenames derived from the reads directory
    :return: suffix of the assembly
    :rtype: string
    :raises FileNotFoundError: if the corresponding reads cannot be found or are not complete

    """

    assembly_suf = ''

    for fa in os.listdir(assembly_dir):
        basename, assembly_suf = fa.rsplit('.', 1)
        if basename not in basenames:
            logging.exception('Cannot find the corresponding reads for {}'.format(basename))
            raise FileNotFoundError('Cannot find the corresponding reads for {}'.format(basename))

    return assembly_suf


def check_reads(reads_dir):
    """
    Check if the assembly provided can match corresponding reads pairs

    :param str reads_dir: directory containing metagenomic paired-end reads
    :return: basenames of input files, suffix of the reads
    :rtype: tuple
    :raises FileNotFoundError: if the reads provided in the input directory are not in pairs

    """
    basenames = []
    reads_suf = ''
    interleaved_reads = None
    pattern = re.compile(r'(.*)(\.fq|\.fastq)(\.gz)?')
    for read in os.listdir(reads_dir):
        match_res = pattern.match(read)
        reads_suf = ''.join(match_res.groups('')[1:])

        if match_res.group(1).endswith('_1'):
            if interleaved_reads is True:
                logging.exception('Reads formats are not consistent with each other')
                raise FileNotFoundError('Reads formats are not consistent with each other')
            interleaved_reads = False
            basename = match_res.group(1)[:-2]
            basenames.append(basename)
            if not os.path.exists(os.path.join(reads_dir, basename) + '_2' + reads_suf):
                logging.exception('Failed to find the corresponding reverse reads to {}'.format(read))
                raise FileNotFoundError('Failed to find the corresponding reverse reads to {}'.format(read))
        elif match_res.group(1).endswith('_2'):
            pass
        else:
            if interleaved_reads is False:
                logging.exception('Reads formats are not consistent with each other')
                raise FileNotFoundError('Reads formats are not consistent with each other')
            interleaved_reads = True
            basename = match_res.group(1)
            basenames.append(basename)

    return basenames, reads_suf, interleaved_reads
