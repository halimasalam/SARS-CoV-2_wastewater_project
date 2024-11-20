# This pipeline requires installing artic bioinformatic tool [https://github.com/artic-network/fieldbioinformatics/] and runs in the artic environment

# Also, the sample sheet is a critical component of the pipeline, as it serves as a guide for determining which samples should undergo analysis. The most vital fields are: "sample_id" and "sequencing_sample_id"

# import required libraries
import os
import subprocess
import argparse
import pandas as pd

# Function to install required packages
def install_packages():
    print(f"Installing packages")
    packages = ["pigz", "fastp", "ivar", "varscan"] 
    for package in packages:
        command = f"mamba install -y -c bioconda {package}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE)
        print(f"Installed package: {package}")

# Function to read the sample sheet and extract sample information
def read_sample_sheet(filename):
    samples = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split(',')
        for line in lines[1:]:
            fields = line.strip().split(',')
            sample = dict(zip(header, fields))
            samples.append(sample)
    return samples

# Function to run pigz compression on Illumina reads
def run_pigz(sample, readsdir_ill, outdir, split_cpus):
    sequencing_sample_id = sample['sequencing_sample_id']
    print(f"Running pigz for sample {sequencing_sample_id}")
    # Run Pigz on the R1 and R2 Illumina reads and save them in the output directory
    subprocess.run(f"pigz -{split_cpus} -cd {readsdir_ill}/{sequencing_sample_id}-R1-001.fastq.gz | pigz -{split_cpus} > {outdir}/pigz/{sequencing_sample_id}.collapsed-R1.fastq.gz", shell=True, stdout=subprocess.PIPE)
    subprocess.run(f"pigz -{split_cpus} -cd {readsdir_ill}/{sequencing_sample_id}-R2-001.fastq.gz | pigz -{split_cpus} > {outdir}/pigz/{sequencing_sample_id}.collapsed-R2.fastq.gz", shell=True, stdout=subprocess.PIPE)

# Function to run Fastp for quality trimming of Illumina reads
def run_fastp(sample, outdir, split_cpus):
    sequencing_sample_id = sample['sequencing_sample_id']
    print(f"Running fastp for sample {sequencing_sample_id}") 
    # Run fastp to trim and filter Illumina reads
    subprocess.run(f"fastp \
    -i {outdir}/pigz/{sequencing_sample_id}.collapsed-R1.fastq.gz \
    -I {outdir}/pigz/{sequencing_sample_id}.collapsed-R2.fastq.gz \
    -o {outdir}/fastp/{sequencing_sample_id}.trimmed-R1.fastq.gz \
    -O {outdir}/fastp/{sequencing_sample_id}.trimmed-R2.fastq.gz \
    -t {split_cpus} \
    -j {outdir}/fastp/{sequencing_sample_id}.json", shell=True, stdout=subprocess.PIPE) 

# Function to run BWA-MEM alignment for Illumina reads
def run_bwa_mem(sample, outdir, split_cpus, fasta):
    sequencing_sample_id = sample['sequencing_sample_id']
    print(f"Running bwa_mem for sample {sequencing_sample_id}")
    # Index the reference fasta and align the trimmed Illumina reads
    subprocess.run(f"bwa index {fasta}", shell=True, stdout=subprocess.PIPE)
    subprocess.run(f"bwa mem \
    -t {split_cpus} \
    {fasta} \
    {outdir}/fastp/{sequencing_sample_id}.trimmed-R1.fastq.gz \
    {outdir}/fastp/{sequencing_sample_id}.trimmed-R2.fastq.gz \
    | samtools view -@ {split_cpus} -bh -o {outdir}/bwa_mem/{sequencing_sample_id}.bam", shell=True, stdout=subprocess.PIPE)

# Function to run ivar variant calling for Illumina reads
def run_ivar(sample, outdir, nimagen_scheme):
    sequencing_sample_id = sample['sequencing_sample_id']
    print(f"Running ivar for sample {sequencing_sample_id}")
    # Sort, index, and trim Illumina BAM file, then run ivar
    command = (
        f"samtools sort -o {outdir}/bwa_mem/{sequencing_sample_id}.sorted.bam {outdir}/bwa_mem/{sequencing_sample_id}.bam && "
        f"samtools index {outdir}/bwa_mem/{sequencing_sample_id}.sorted.bam && "
        f"ivar trim -i {outdir}/bwa_mem/{sequencing_sample_id}.sorted.bam -b {nimagen_scheme} -p {outdir}/ivar/{sequencing_sample_id}.trimmed.bam && "
        f"samtools sort -o {outdir}/ivar/{sequencing_sample_id}.sorted.trimmed.bam {outdir}/ivar/{sequencing_sample_id}.trimmed.bam && "
        f"samtools index {outdir}/ivar/{sequencing_sample_id}.sorted.trimmed.bam")
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)

# Function to create a modified primer scheme for artic pipeline
def artic_nimagen_primer_scheme(artic_scheme_folder, nimagen_scheme):
    print(f"Create artic_nimagen_primer_scheme")
    subprocess.run(f"cp -R {artic_scheme_folder}/nCoV-2019/V3 {artic_scheme_folder}/nCoV-2019/V_nimagen", shell=True, stdout=subprocess.PIPE)
    # Modify and create a new primer scheme file for artic pipeline
    command = f"awk -F'\t' 'BEGIN {{OFS=\"\t\"}} {{gsub(/right/, \"RIGHT\", $4); gsub(/left/, \"LEFT\", $4); if ($5 == 1) $5 = \"nCoV-2019_1\"; if ($5 == 2) $5 = \"nCoV-2019_2\"; NF--}}1' {nimagen_scheme} > {artic_scheme_folder}/nCoV-2019/V_nimagen/nCoV-2019.scheme.bed" 
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

# Function to merge ONT fastq files if pass and fail folders exist
def merge_ont_fastq(ont_dir_path, ont_dir):
    pass_dir = os.path.join(ont_dir_path, "fastq_pass")
    fail_dir = os.path.join(ont_dir_path, "fastq_fail")
     # Merge ONT fastq files if both pass and fail folders exist
    if os.path.exists(pass_dir) and os.path.exists(fail_dir):
        merged_fastq_path = os.path.join(ont_dir_path, f"catted_{ont_dir}.fastq")
        subprocess.run(f'cat "{pass_dir}"/*.fastq.gz "{fail_dir}"/*.fastq.gz > "{merged_fastq_path}"', shell=True)
        print(f"Fastq files merged for {ont_dir}")
    else:
        print(f"Skipping merge for {ont_dir}: Missing pass and/or fail folders")

# Function to run artic guppyplex for ONT reads
def artic_guppyplex(ont_dir_path, ont_dir, split_cpus):
    print(f"Running artic_guppyplex for sample {ont_dir}")
    # Run artic guppyplex to prepare ONT reads for artic minion
    subprocess.run(f"artic guppyplex \
    --min-length 80 \
    --max-length 480 \
    --directory {ont_dir_path} \
    --prefix guppyplex", shell=True, stdout=subprocess.PIPE)

# Function to run artic minion for ONT reads
def artic_minion(artic_scheme_folder, ont_dir, ont_dir_path, split_cpus):
    print(f"Running artic_guppyplex for sample {ont_dir}")
    # Run artic minion using ONT reads and specific scheme
    subprocess.run(f"artic minion \
    --normalise 1000000 \
    --threads {split_cpus} \
    --medaka --medaka-model r1041_e82_400bps_sup_g615 \
    --scheme-directory {artic_scheme_folder} \
    --read-file {ont_dir_path}/guppyplex_{ont_dir}.fastq \
    nCoV-2019/V_nimagen {ont_dir}", shell=True, stdout=subprocess.PIPE)

# Function to match ONT nimagen indexes with corresponding illumina well plates
def match_nimagen_indexes_with_well_plates(nimagen_index, ont_dir, ont_directories):
    csv = pd.read_csv(nimagen_index, sep="\t")
    csv = csv.rename(columns={"Index combination name": "Index_combination_name"})
    csv = csv.rename(columns={"Plate location": "Plate_location"})
    csv["index_number"] = csv.Index_combination_name.str[-4:7]
    nimagen_plates = {dir[-4:]: row['Plate_location'] for dir in ont_directories for index, row in csv.iterrows() if row['index_number'] == dir[-4:]}
    Nimagen_plates = {v: k for k, v in nimagen_plates.items()}
    return Nimagen_plates

# Function to create folders for analysis result with naming format of [ONT nimagen indexes with corresponding illumina well plates]
def create_analysis_directory(nimagen_index, ont_dir, outdir, ont_directories):
    print("Matching nimagen indexes with well-plates")
    Nimagen_plates = match_nimagen_indexes_with_well_plates(nimagen_index, ont_dir, ont_directories)
    illumina_bams = [f for f in os.listdir(f"{outdir}/ivar/") if f.endswith(".bam")]
    
    for bam in illumina_bams:
        plate = Nimagen_plates.get(bam.split('-')[0])
        if plate:
            os.makedirs(f"{outdir}/analysis/{bam.split('.')[0]}-{plate[-4:]}", exist_ok=True)

# Function to run samtools mpileup and mosdepth for illumina reads
def samtools_and_mosdepth_ill(outdir, fasta):
    for bam in [f for f in os.listdir(f"{outdir}/ivar/") if f.endswith(".sorted.trimmed.bam")]:
        for joint_folder in [f for f in os.listdir(f"{outdir}/analysis/") if bam.split('-')[0] == f.split('-')[0]]:
            joint_folder_path = os.path.join(f"{outdir}/analysis/", joint_folder)
            print(f"Running samtools mpileup for {joint_folder} Illumina")
            subprocess.run(f"samtools mpileup \
            -f {fasta} \
            -q 10 \
            -d 1000000 \
            {outdir}/ivar/{bam} > {joint_folder_path}/{joint_folder}.ill.pileup",
            shell=True, stdout=subprocess.PIPE)

            print(f"Running mosdepth for {joint_folder} illumina")
            subprocess.run(f"mosdepth \
            --no-per-base  {joint_folder_path}/{joint_folder}.ill {outdir}/ivar/{bam}",
            shell=True, stdout=subprocess.PIPE)

            print(f"Running samtools depth for {joint_folder} illumina")
            subprocess.run(f"samtools depth \
            -aa {outdir}/ivar/{bam} >  {joint_folder_path}/{joint_folder}.ill.samtools.cov",
            shell=True, stdout=subprocess.PIPE)

# Function to run samtools mpileup and mosdepth for ont reads
def samtools_and_mosdepth_ont(outdir, ont_dir_path, fasta):  
    for bam in [f for f in os.listdir(ont_dir_path) if f.endswith(".primertrimmed.rg.sorted.bam")]:
        for joint_folder in [f for f in os.listdir(f"{outdir}/analysis/") if bam.split(".")[0][-4:] == f[-4:]]:
            joint_folder_path = os.path.join(f"{outdir}/analysis/", joint_folder)
            print(f"Running samtools mpileup for {joint_folder} ONT")
            subprocess.run(f"samtools mpileup \
            -f {fasta} \
            -q 10 \
            -d 1000000 \
            {ont_dir_path}/{bam} > {joint_folder_path}/{joint_folder}.ont.pileup",
            shell=True, stdout=subprocess.PIPE)
    
            print(f"Running mosdepth for {joint_folder} ONT")
            subprocess.run(f"mosdepth \
            --no-per-base  {joint_folder_path}/{joint_folder}.ont {ont_dir_path}/{bam}",
            shell=True, stdout=subprocess.PIPE)
    
            print(f"Running samtools depth for {joint_folder} ONT")
            subprocess.run(f"samtools depth \
            -aa {ont_dir_path}/{bam} >  {joint_folder_path}/{joint_folder}.ont.samtools.cov",
            shell=True, stdout=subprocess.PIPE)

# Function to varscan for illumina reads
def varscan_baseline_and_optimal_ill(outdir):
    for joint_folder in [f for f in os.listdir(f"{outdir}/analysis/")]:
        joint_folder_path = os.path.join(f"{outdir}/analysis/", joint_folder)
        for pileup in [f for f in os.listdir(joint_folder_path) if f.endswith(".ill.pileup")]:
            print(f"Running varscan at baseline threshold for {joint_folder} illumina")
            subprocess.run(f"varscan mpileup2snp \
            {joint_folder_path}/{pileup} \
            --min-var-freq 0.01 \
            --p-value 1 \
            --min-coverage 1 \
            --strand-filter 0 \
            --min-reads2 1 > {joint_folder_path}/{joint_folder}.ill.varscan.baseline.tsv",
            shell=True, stdout=subprocess.PIPE)

            print(f"Running varscan at optimal threshold for {joint_folder} illumina")
            subprocess.run(f"varscan mpileup2snp \
            {joint_folder_path}/{pileup} \
            --min-var-freq 0.06 \
            --min-coverage 4 \
            --strand-filter 0 > {joint_folder_path}/{joint_folder}.ill.varscan.optimal.tsv",
            shell=True, stdout=subprocess.PIPE)

# Function to run varscan for ONT reads
def varscan_baseline_and_optimal_ont(outdir):
    for joint_folder in [f for f in os.listdir(f"{outdir}/analysis/")]:
        joint_folder_path = os.path.join(f"{outdir}/analysis/", joint_folder)
        for pileup in [f for f in os.listdir(joint_folder_path) if f.endswith(".ont.pileup")]:
            print(f"Running varscan at baseline threshold for {joint_folder} ONT with no strand bias ")
            subprocess.run(f"varscan mpileup2snp \
            {joint_folder_path}/{pileup} \
            --min-var-freq 0.01 \
            --p-value 1 \
            --min-coverage 1 \
            --strand-filter 0 \
            --min-reads2 1 > {joint_folder_path}/{joint_folder}.ont.varscan.nosb.baseline.tsv",
            shell=True, stdout=subprocess.PIPE)

            print(f"Running varscan at baseline threshold for {joint_folder} ONT with default strand bias filter")
            subprocess.run(f"varscan mpileup2snp \
            {joint_folder_path}/{pileup} \
            --min-var-freq 0.01 \
            --p-value 1 \
            --min-coverage 1 \
            --min-reads2 1 > {joint_folder_path}/{joint_folder}.ont.varscan.baseline.tsv",
            shell=True, stdout=subprocess.PIPE)

            print(f"Running varscan at optimal threshold for {joint_folder} ONT")
            subprocess.run(f"varscan mpileup2snp \
            {joint_folder_path}/{pileup} \
            --min-var-freq 0.06 \
            --min-coverage 4 > {joint_folder_path}/{joint_folder}.ont.varscan.optimal.tsv",
            shell=True, stdout=subprocess.PIPE)
            
# Main function to manage the entire workflow
def main(args):
    install_packages()
    # Create necessary output directories
    os.makedirs(f"{args.outdir}/pigz", exist_ok=True)
    os.makedirs(f"{args.outdir}/fastp", exist_ok=True)
    os.makedirs(f"{args.outdir}/bwa_mem", exist_ok=True)
    os.makedirs(f"{args.artic_scheme_folder}/nCoV-2019/V_nimagen", exist_ok=True)
    os.makedirs(f"{args.outdir}/ivar", exist_ok=True)

    # Read sample sheet to get sample information
    samples = read_sample_sheet(args.sample_sheet)
    fasta = f"{args.artic_scheme_folder}/nCoV-2019/V_nimagen/nCoV-2019.reference.fasta"
    
    for sample in samples:
        run_pigz(sample, args.readsdir_ill, args.outdir, args.split_cpus)
    
    for sample in samples:
        run_fastp(sample, args.outdir, args.split_cpus)
    
    for sample in samples:
        run_bwa_mem(sample, args.outdir, args.split_cpus, fasta)
    
    for sample in samples:
        run_ivar(sample, args.outdir, args.nimagen_scheme)

    artic_nimagen_primer_scheme(args.artic_scheme_folder, args.nimagen_scheme)
    
    ont_directories = [d for d in os.listdir(args.readsdir_ont)]
    for ont_dir in ont_directories:
        ont_dir_path = os.path.join(args.readsdir_ont, ont_dir)
        os.chdir(ont_dir_path)
        merge_ont_fastq(args.readsdir_ont, ont_dir)
        artic_guppyplex(ont_dir_path, ont_dir, args.split_cpus)
        artic_minion(args.artic_scheme_folder, ont_dir, ont_dir_path,args.split_cpus)

    create_analysis_directory(args.nimagen_index, ont_dir,args.outdir, ont_directories)

    samtools_and_mosdepth_ill(args.outdir, fasta)
    
    for ont_dir in ont_directories:
        ont_dir_path = os.path.join(args.readsdir_ont, ont_dir)
        samtools_and_mosdepth_ont(args.outdir, ont_dir_path, fasta)

    varscan_baseline_and_optimal_ill(args.outdir)

    varscan_baseline_and_optimal_ont(args.outdir)

 # Check if the script is being run directly       
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process Illumina and ONT fastq files.")
    parser.add_argument("--readsdir_ill", required=True, help="Path to raw illumina fastq files directory")
    parser.add_argument("--readsdir_ont", required=True, help="Path to raw ont fastq files directory")
    parser.add_argument("--artic_scheme_folder", required=True, help="Path to artic_primer_scheme folder")
    parser.add_argument("--nimagen_scheme", required=True, help="Path to nimagen_primer_scheme")
    parser.add_argument("--nimagen_index", required=True, help="Path to Nimagen IDX index file ")
    parser.add_argument("--sample-sheet", required=True, help="Path to sample sheet file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--split-cpus", type=int, default=4, help="Number of split CPUs")
    args = parser.parse_args()
    
    # Call the main function
    main(args)
