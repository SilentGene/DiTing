<img src="./DiTing_logo.jpg" width="250" height="250">

# DiTing
## Etymology
**DiTing** is a Chinese mythical creature who knows everything when he puts ears on the earth's surface. Parallelly, this program is developed to recognize biogeochemical cycles from environmental omic data accurately and efficiently.    
**谛听(DiTing)** 若伏在地下，一霎时，便可将四大部洲山川社稷、洞天福地之间， 蠃虫、鳞虫、毛虫、羽虫、昆虫，天仙、地仙、神仙、人仙、鬼仙，顾鉴善恶，察听贤愚。

## Citation
To cite DiTing please use  
> Xue CX, Lin H, Zhu XY, Liu J, Zhang Y, Rowley G, Todd JD, Li M, Zhang XH. DiTing: A Pipeline to Infer and Compare Biogeochemical Pathways From Metagenomic and Metatranscriptomic Data. Front Microbiol. 2021 Aug 2;12:698286. doi: 10.3389/fmicb.2021.698286.    

## Introduction
DiTing is designed to determine the relative abundance of metabolic and biogeochemical functional pathways in a set of given metagenomic/metatranscriptomic data. The input is expected to be a folder containing a group of paired-end clean reads. These reads will be assembled, annotated, and parsed for producing a table of relative abundance of elemental/biogeochemical cycling pathways (e.g., Nitrogen, Carbon, Sulfur) in each sample. Sketch maps and heatmaps will also be produced accordingly for comparing biogeochemical functions visually.

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

### Option 1: Conda Environment Setup (Recommended)
You can build the completely reproducible stack straight from the `environment.yaml`:

```bash
# 1. Download the repo
git clone https://github.com/xuechunxu/DiTing.git
cd DiTing

# 2. Build the exact environment mapping resolving all snakemake dependencies
conda env create -f environment.yaml

# 3. Activate the environment
conda activate diting-env
```
Once inside the environment, the `diting` CLI command will be available natively as an entry point. 

### Option 2: Build from Source
If setting up via pure pip, you need to ensure the `Dependencies` outlined above are already correctly placed onto your system PATH, then:
```bash
git clone https://github.com/xuechunxu/DiTing.git
cd DiTing
pip install -e .
```

#### Database Downloads
DiTing requires [KofamKOALA hmm database](https://www.genome.jp/tools/kofamkoala/). You can download and extract the database using the following command:

```bash
diting-download-db -o <kofam_database>
```
This will download `ko_list` and the profile HMMs into the specified directory.

## Running
### 1. One step running
Instead of invoking `python diting.py`, the CLI automatically handles forming connections into the `Snakemake` pipeline.

```bash
# 1. Download database
diting-download-db -o kofam_database

# 2. Run dieting
diting -r <clean_reads_dir> -o <output_dir> -p kofam_database/profiles -k kofam_database/ko_list
diting -r <clean_reads_dir> -a <metagenomic_assembly> -o <output_dir> -p kofam_database/profiles -k kofam_database/ko_list
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
The input is the `<clean_reads_dir>` folder containing a group of paired-end metagenomic clean reads, looks like: 
```
sample_one_1.fastq
sample_one_2.fastq
...
sample_three_1.fastq
sample_three_2.fastq
```
The paired-end metagenomic clean reads should end with `.fq`, `.fq.gz`, `.fastq`, or `.fastq.gz`.
The interleaved reads are also supported.

### 2. Optional parameter
#### 2.1 --spades
Using `metaSPAdes` instead of `megahit` to assemble reads

Consider setting memory limitation by `-m` when usign `SPAdes` as assembler

`-m(--memory) <int>` default: 50 (in Gb)

#### 2.2 -a (--assembly) metagenomic assembly
Path to a folder containing metagenomic assemblies corresponding to the provided reads, which is expected to have the same base name as the reads. The reads will not be assembled when this parameter was used.

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
The intermediate `.sam` files would be retained if this flag was used. 
```bash
diting -r <clean_reads_dir> -o <output_dir> -n 12 --noclean -p <profiles_dir> -k <ko_list>
```
#### 2.6 -vis (--visualization) pathways_relative_abundance.tab
Visualization can also be executed independently, which allows users to adjust the final result table (e.g., merge some similar samples) before the visualization.
```bash
diting -vis <pathways_relative_abundance.tab>
```

#### 2.7 --dry-run
Perform a dry run of the snakemake DAG pipeline to view exactly which sequence tasks execute.
```bash
diting -r <clean_reads_dir> -o <output_dir> --dry-run -p <profiles_dir> -k <ko_list>
```
### 3. Output
#### 3.1 Table
- `pathways_relative_abundance.tab` :The final result with the relative abundance of pathways in each sample. 
- `ko_abundance_among_samples.tab` : A table with the relative abundance of each `k_number` of KEGG annotation is produced in `KEGG_annotation` folder. 

#### 3.2 Visualization
- `carbon_cycle_sketch.png`, `nitrogen_cycle_sketch.png`, `DMSP_cycle_sketch.png` and `sulfur_cycle_sketch.png`
Sketch maps regarding carbon, nitrogen and sulfur cycles
- `carbon_cycle_heatmap.pdf`, `nitrogen_cycle_heatmap.pdf`, `sulfur_cycle_heatmap.pdf` and `other_cycle_heatmap.pdf`
Heatmaps regarding carbon, nitrogen, sulfur cycles and other pathways

Example:
`sketch`look like:
<img src="./example/diting.out/sketch.png" width="792" height="624">

`heatmap`look like:
<img src="./example/diting.out/heatmap.png" width="792" height="627">

## Copyright
Xue Chunxu, xuechunxu (at) outlook.com  
Heyu Lin, heyu.lin (at) qut.edu.au  
Xiaoyu Zhu, xiaoyuzhu321 (at) 126.com  
Xiao-Hua Zhang, xhzhang (at) ouc.edu.cn  
Lab of Microbial Oceanography  
College of Marine Life Sciences, Ocean University of China, Qingdao 266003, China
