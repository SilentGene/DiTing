import argparse
import os
import sys
import urllib.request
import gzip
import tarfile
import shutil
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser(description="Download and extract KEGG kofam database for DiTing")
    parser.add_argument('-o', '--outdir', metavar='output_dir', dest='o', type=str, required=True, help='Output directory for the database')
    
    args = parser.parse_args()
    
    out_dir = os.path.abspath(args.o)
    if not os.path.exists(out_dir):
        logging.info(f"Creating directory: {out_dir}")
        os.makedirs(out_dir, exist_ok=True)
        
    url_ko_list = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
    url_profiles = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
    
    path_ko_list_gz = os.path.join(out_dir, 'ko_list.gz')
    path_ko_list = os.path.join(out_dir, 'ko_list')
    path_profiles_tar_gz = os.path.join(out_dir, 'profiles.tar.gz')
    
    try:
        logging.info("Downloading ko_list.gz...")
        with urllib.request.urlopen(url_ko_list) as response, open(path_ko_list_gz, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
            
        logging.info("Downloading profiles.tar.gz...")
        with urllib.request.urlopen(url_profiles) as response, open(path_profiles_tar_gz, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
            
        logging.info("Decompressing ko_list.gz...")
        with gzip.open(path_ko_list_gz, 'rb') as f_in, open(path_ko_list, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
        logging.info("Extracting profiles.tar.gz...")
        with tarfile.open(path_profiles_tar_gz, "r:gz") as tar:
            tar.extractall(path=out_dir)
            
        # Clean up
        os.remove(path_ko_list_gz)
        os.remove(path_profiles_tar_gz)
        
        logging.info(f"Database has been downloaded and extracted successfully to {out_dir}")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
