import os
import subprocess

input_bed12 = snakemake.input.bed12
input_gtf = snakemake.input.gtf
input_reads = snakemake.input.reads
output_gtf = snakemake.output.gtf
params_reference_fa = snakemake.input.ref_fa
params_output_prefix = snakemake.params.output_prefix
params_bed_split_size = snakemake.params.bed_split_size
log = snakemake.log[0]
threads = snakemake.threads
tissue = snakemake.wildcards.tissue

with open(log, "w") as log_file:
    if os.path.getsize(input_bed12) > float(params_bed_split_size):
        base_dir = os.path.dirname(input_bed12)
        # create {base_dir}/split_{tissue}/chr{chrom}.bed files
        split_dir = os.path.join(base_dir, f"split_{tissue}")
        os.makedirs(split_dir, exist_ok=True)

        # Use awk to split the BED file by chromosome
        awk_command = f"awk '{{print > \"{split_dir}/\"$1\".bed\"}}' {input_bed12}"
        subprocess.run(awk_command, shell=True, check=True)

        output_dir = os.path.join(params_output_prefix, f"{tissue}_split")
        os.makedirs(output_dir, exist_ok=True)
        # Run FLAIR on each chromosome
        for chrom_bed in os.listdir(split_dir):
            chrom = chrom_bed.split(".")[0]
            output_prefix = os.path.join(output_dir, chrom)
            output_gtf = f"{output_prefix}.isoforms.gtf"
            if os.path.exists(output_gtf):
                log_file.write(f"Skipping split FLAIR collapse on {chrom} as {output_gtf} already exists\n")
                continue
            log_file.write(f"Running split FLAIR collapse on {chrom}\n")
            subprocess.run(f"flair collapse --query {os.path.join(split_dir, chrom_bed)} \
                           --genome {params_reference_fa} --reads {input_reads} --gtf {input_gtf} \
                           --threads {threads} --output {output_prefix}", shell=True, check=True,
                           stdout=log_file, stderr = subprocess.STDOUT, text=True)

        log_file.write(f"Concatenating split files\n")
        log_file.write(f"Found files: {os.listdir(split_dir)}\n")
        # Concatenate all the collapsed files
        collapsed_files = " ".join([os.path.join(output_dir, f"{chrom}.isoforms.gtf") for chrom in os.listdir(split_dir)])
        subprocess.run(f"cat {collapsed_files} > {output_gtf}", shell=True, check=True,
                       stdout=log_file, stderr = subprocess.STDOUT, text=True)
    else:
        log_file.write(f"Running FLAIR collapse on {input_bed12}\n")
        subprocess.run(f"flair collapse --query {input_bed12} --genome {params_reference_fa} \
                       --reads {input_reads} --gtf {input_gtf} --threads {threads} \
                       --output {params_output_prefix}", shell=True, check=True, stdout=log_file,
                       stderr = subprocess.STDOUT, text=True)
