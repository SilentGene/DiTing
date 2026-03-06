import argparse
import os
import sys
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="DiTing: A Pipeline to Infer and Compare Biogeochemical Pathways")
    parser.add_argument('-r', '--reads', metavar='input_reads', dest='r', type=str, required=False, help='folder containing reads to be used as input')
    parser.add_argument('-o', '--outdir', metavar='output_dir', dest='o', type=str, required=False, help='output directory')
    parser.add_argument('-a', '--assembly', metavar='metagenomic_assembly', dest='a', type=str, help='folder containing metagenomic assemblies corresponding to provided reads, which should have the same basename as the reads')
    parser.add_argument('-n', '--threads', metavar='threads', dest='n', type=int, default=4, help='threads that will be used')
    parser.add_argument('--noclean', dest='nc', action='store_true', default=False, help='The sam files would be retained if this flag is used')
    parser.add_argument('-vis', '--visualization', metavar='pathways_relative_abundance.tab', dest='vis', type=str, default=False, help='A table for visualization')
    parser.add_argument('--spades', dest='spades', action='store_true', default=False, help='metaSPAdes will be used for assembling instead of megahit if this flag is used')
    parser.add_argument('-m', '--memory', metavar='memory', dest='m', type=int, default=50, help='Memory that will be used by metaSPAdes (in Gb). Default=50G')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run of the snakemake pipeline')
    parser.add_argument('--snakemake-args', nargs=argparse.REMAINDER, help='Additional arguments to pass to snakemake')
    
    args = parser.parse_args()

    # Find the Snakefile
    snakefile = Path(__file__).resolve().parent / "workflow" / "Snakefile"
    if not snakefile.exists():
        print(f"Error: Snakefile not found at {snakefile}", file=sys.stderr)
        sys.exit(1)

    cmd = [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile), "--cores", str(args.n)]
    
    if args.dry_run:
        cmd.append("-n")
    
    # Pass configuration parameters
    config_args = []
    if args.vis:
        config_args.append(f"vis={os.path.abspath(args.vis)}")
    else:
        if not args.r:
            print("Error: The following arguments are required: -r/--reads", file=sys.stderr)
            sys.exit(1)
        if not args.o:
            print("Error: The following arguments are required: -o/--outdir", file=sys.stderr)
            sys.exit(1)
            
        config_args.extend([
            f"reads_dir={os.path.abspath(args.r)}",
            f"out_dir={os.path.abspath(args.o)}",
            f"spades={args.spades}",
            f"memory={args.m}",
            f"threads={args.n}",
            f"noclean={args.nc}"
        ])
        if args.a:
            config_args.append(f"assembly_dir={os.path.abspath(args.a)}")
            
    if config_args:
        cmd.append("--config")
        cmd.extend(config_args)
        
    if args.snakemake_args:
        cmd.extend(args.snakemake_args)

    print(f"Running Snakemake command: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
