<img src="./DiTing_logo.jpg" width="250" height="250">

# DiTing
![Version](https://img.shields.io/badge/version-v2.0.1-blue)
![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)
![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.0-brightgreen)

> [!NOTE]
> **🚀 Major Update (v2.0):** DiTing has been boosted to version 2! Compared to the initially published version, this release introduces several significant upgrades:
> 1. **Complete Architecture Rewrite**: The entire pipeline has been rewritten using `Snakemake`, providing robust workflow management, better parallelization, and the ability to resume execution from breakpoints.
> 2. **Enhanced Annotation Engine**: We have replaced the manual `hmmsearch` parsing system with the standardized `kofamscan` engine for KEGG annotations. This effectively resolves the frequent parsing errors and compatibility issues previously encountered with raw `hmmsearch` outputs.
> 3. **Maintenance Notice**: DiTing is no longer under highly active development. However, if you encounter any issues during usage, please feel free to open an issue on GitHub, and we will do our best to address it!

## Etymology
**DiTing** is a Chinese mythical creature who knows everything when he puts his ears to the earth's surface. Similarly, this program is developed to accurately and efficiently recognize biogeochemical cycles from environmental omic data.    
**谛听 (DiTing)** 若伏在地下，一霎时，便可将四大部洲山川社稷、洞天福地之间，蠃虫、鳞虫、毛虫、羽虫、昆虫，天仙、地仙、神仙、人仙、鬼仙，顾鉴善恶，察听贤愚。

## Citation
To cite DiTing please use  
> Xue CX, Lin H, Zhu XY, Liu J, Zhang Y, Rowley G, Todd JD, Li M, Zhang XH. DiTing: A Pipeline to Infer and Compare Biogeochemical Pathways From Metagenomic and Metatranscriptomic Data. Front Microbiol. 2021 Aug 2;12:698286. doi: [10.3389/fmicb.2021.698286](https://doi.org/10.3389/fmicb.2021.698286).    

## Introduction
DiTing is designed to determine the relative abundance of metabolic and biogeochemical functional pathways in a set of given metagenomic or metatranscriptomic data. The input should be a folder containing a group of paired-end clean reads. These reads will be assembled, annotated, and parsed to produce a table detailing the relative abundance of elemental and biogeochemical cycling pathways (e.g., Nitrogen, Carbon, Sulfur, and DMSP) in each sample. Sketch maps and heatmaps will also be produced to visually compare these biogeochemical functions.

## Procedure
![image](./Flow_chart.png)

## Dependencies
DiTing now relies on **Snakemake >= 9.0** to orchestrate the pipeline and uses `hatchling` to package the modules. The underlying bioinformatics dependencies remain:
* [Megahit](https://github.com/voutcn/megahit)
* [SPAdes](https://cab.spbu.ru/software/spades/)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [bwa](https://github.com/lh3/bwa)
* [BBMap](https://github.com/BioInfoTools/BBMap)
* [HMMER3](http://hmmer.org/)
* [KofamScan](https://github.com/takaram/kofam_scan)
* Python modules (handled strictly by conda/pip): 
    * `pandas`, `matplotlib`, `opencv-python`, `Pillow`, `seaborn`
* KofamKOALA hmm database (ftp://ftp.genome.jp/pub/db/kofam/)

## Installation
Recommended configuration:  
```
CPU threads ≥ 8  
RAM ≥ 64 Gb
```

### Conda Environment Setup

```bash
# 1. Download the repo
git clone https://github.com/SilentGene/DiTing.git
cd DiTing

# 2. Build the conda environment called 'diting'
conda env create -f environment.yaml

# 3. Activate the environment
conda activate diting
```

#### Database Downloads
DiTing requires [KofamKOALA hmm database](https://www.genome.jp/tools/kofamkoala/). You can download and extract the database using the following command:

```bash
diting-download-db -o <kofam_database>
```
This will download `ko_list` and the profile HMMs folder `profiles` into the specified directory. Alternatively, you can download the database manually and extract it into the specified directory.

```bash
mkdir kofam_database
cd kofam_database
wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz 
wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz 
gzip -d ko_list.gz
tar zxvf profiles.tar.gz 
```

## Running
### 1. One step running

```bash
# from reads (interleaved or paired-end)
diting -r <clean_reads_dir> -o <output_dir> -p kofam_database/profiles -k kofam_database/ko_list

# from reads and assembly (contigs)
diting -r <clean_reads_dir> -a <metagenomic_assembly_dir> -o <output_dir> -p kofam_database/profiles -k kofam_database/ko_list
```
Example reads run:  
```bash
#download the example reads  
Google Drive:  
URL: https://drive.google.com/file/d/132605rtKuA-Xx--eh3aC7i5WIExNWl5k/view?usp=sharing
after download, run:
unzip Clean-reads_interleaved.zip

# run Example
diting -r Clean-reads_interleaved -o Clean-reads_interleaved.diting.out -p kofam_database/profiles -k kofam_database/ko_list
```
The input is the `<clean_reads_dir>` folder containing a group of paired-end metagenomic clean reads, formatted as follows:
```
sample_one_1.fastq
sample_one_2.fastq
...
sample_three_1.fastq
sample_three_2.fastq
```
The paired-end metagenomic clean reads should end with `.fq`, `.fq.gz`, `.fastq`, or `.fastq.gz`.
Interleaved reads are also supported and should be formatted as follows:
```
sample_one.fq.gz
sample_two.fq.gz
sample_three.fq.gz
```

### 2. Optional parameter
#### 2.1 --spades
Using `metaSPAdes` instead of `megahit` to assemble reads

Consider setting memory limitation by `-m` when usign `SPAdes` as assembler

`-m(--memory) <int>` default: 50 (in Gb)

#### 2.2 -a (--assembly) metagenomic assembly
Path to a folder containing metagenomic assemblies corresponding to the provided reads, which is expected to have the same base names as the reads. The reads will not be assembled when this parameter is used.

```bash
diting -r <clean_reads_dir> -a <metagenomic_assembly> -o <output_dir> -p <profiles_dir> -k <ko_list>
```
The `<metagenomic_assembly>` folder looks like: 
```
sample_one.fa
sample_two.fa
sample_three.fa
```

#### 2.4 -n (--threads) number of threads
Number of threads to run (default: 4)

```bash
diting -r <clean_reads_Dir> -a <metagenomic_assembly> -o <output_dir> -n 20 -p <profiles_dir> -k <ko_list>
```
#### 2.5 --noclean
The intermediate `.sam` files will be retained if this flag is used.
```bash
diting -r <clean_reads_dir> -o <output_dir> -n 12 --noclean -p <profiles_dir> -k <ko_list>
```
#### 2.6 -vis (--visualization) pathways_relative_abundance.tab
Visualization can also be executed independently, which allows users to adjust the final result table (e.g., merge some similar samples) before the visualization.
```bash
diting -vis <pathways_relative_abundance.tab>
```

#### 2.7 How to resume a failed run?
As `Snakemake` is used to manage the pipeline, we can resume a failed run by simply running the same command again.
```bash
diting <original_arguments>
```
### 3. Output
#### 3.1 Table
- `pathways_relative_abundance.tab`: The final result containing the relative abundance of pathways in each sample. 
- `ko_abundance_among_samples.tab`: A table containing the relative abundance of each KEGG annotation `k_number`, produced inside the `KEGG_annotation` folder. 

#### 3.2 Visualization
- `carbon_cycle_sketch.png`, `nitrogen_cycle_sketch.png`, `DMSP_cycle_sketch.png` and `sulfur_cycle_sketch.png`
Sketch maps regarding carbon, nitrogen and sulfur cycles
- `carbon_cycle_heatmap.pdf(.png)`, `nitrogen_cycle_heatmap.pdf(.png)`, `sulfur_cycle_heatmap.pdf(.png)` and `other_cycle_heatmap.pdf(.png)`
Heatmaps regarding carbon, nitrogen, sulfur cycles and other pathways

Example:
`sketch` looks like:
<img src="./example/diting.out/sketch.png" width="792" height="624">

`heatmap` looks like:
<img src="./example/diting.out/heatmap.png" width="792" height="627">

## Parameters
```powershell
$ diting --help
usage: diting [-h] [-v] -r input_reads -o output_dir -p profiles_dir -k ko_list [-a metagenomic_assembly] [-n threads] [-m memory]
              [-vis pathways_relative_abundance.tab] [--spades] [--noclean] [--dry-run] [--snakemake-args ...]

DiTing: A Pipeline to Infer and Compare Biogeochemical Pathways

options:
  -h, --help            show this help message and exit
  -v, --version         show program version number and exit
  -r, --reads input_reads
                        folder containing reads to be used as input
  -o, --outdir output_dir
                        output directory
  -p, --profiles profiles_dir
                        folder containing kofam profiles (*.hmm)
  -k, --ko-list ko_list
                        ko_list file
  -a, --assembly metagenomic_assembly
                        folder containing metagenomic assemblies corresponding to provided reads, which should have the same basename as the reads
  -n, --threads threads
                        threads that will be used
  -m, --memory memory   Memory that will be used by metaSPAdes (in Gb). Default=50G
  -vis, --visualization pathways_relative_abundance.tab
                        A table for visualization
  --spades              metaSPAdes will be used for assembling instead of megahit if this flag is used
  --noclean             The sam files would be retained if this flag is used
  --dry-run             Perform a dry run of the snakemake pipeline
  --snakemake-args ...  Additional arguments to pass to snakemake

```

## Contact
Xue Chunxu, xuechunxu (at) outlook.com  
Heyu Lin, heyu.lin (at) qut.edu.au  
Xiaoyu Zhu, xiaoyuzhu321 (at) 126.com  
Xiao-Hua Zhang, xhzhang (at) ouc.edu.cn

Lab of Microbial Oceanography  
College of Marine Life Sciences, Ocean University of China, Qingdao 266003, China

