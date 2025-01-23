[![Snakemake](https://img.shields.io/badge/snakemake-≥8.16.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

This repository contains the code for my master thesis. A few python scripts, to download data files, some notebooks, mostly for experiments, and most importantly a snakemake pipeline that can be used to reproduce the results of the thesis. The pipeline can be run on a local machine or on in cluster environment, but be aware that some of the tools are **very** memory hungry.

# Snakemake Setup

The pipeline was tested using snakemake version 8.16.0 and mamba 1.5.8. Snakemake will install most of the mamba/conda environments automatically, when you run it the first time. However you will need one environment for snakemake itself, from which you run the pipeline, and one for IsoTools, because the latest release isn't available with pip yet.

Snakemake environment:
```bash
mamba create -n snakemake -c bioconda snakemake=8.16.0
```

IsoTools environment:
```bash
mamba create -n isotools pip
conda activate isotools
# Install IsoTools from the repository
git clone git@github.com:HerwigLab/IsoTools2.git
cd IsoTools2
pip install .
```

Additionally you will need to clone SQANTI3 and StringTie 2 to a path of your choice:
```bash
git clone git@github.com:gpertea/stringtie.git
git clone git@github.com:ConesaLab/SQANTI3.git
```
StringTie requires manual building with `make` (with a C++ compiler which supports the C++ 11 standard (GCC 4.8 or newer)):
```bash
cd stringtie
make -j4 release
```

## Configuration

> [!NOTE]
> The configuration is automatically verified at the start of the execution. If you get an error, check the error message for the missing or wrong configuration.

- Open the snakemake [config](./snakemake-pipeline/config/config.yaml).
- Configure the input files
  - paths to the reference annotation and fasta files
  - read data in the form of a tsv sample table
    - table needs to have the columns `sample ID`, `group` and `file`. Samples from different groups are analyzed separately. File is the path to the sample specific file. It can either be a gzipped fastq file or a bam file. The pipeline will automatically detect the file type and adjust.
  - short read data in the form of a tsv sample table, if no short read data is available use `null`
    - table needs to have the columns `sample ID` and `file`. File is the path to the sample specific gzipped fastq reads file.
  - CAGE data again as a tsv table
    - table needs to have the columns `group` and `file`. File is the path to the sample specific gzipped bed file.
  - PolyA site peaks as the path to the gzipped bed file. Files can be downloaded from the [PolyASite 2 database](https://polyasite.unibas.ch/atlas#2).
- Specify which tools you want to use for the analysis. The pipeline supports IsoTools, FLAIR, IsoQuant and StringTie.
  - IsoTools is special here. You can run it multiple times in parallel with different settings. Any tool that starts with `isotools` is detected. The version specific settings are defined below.
- Configure the tools
  - for FLAIR you can specify additional command line parameters for their `correct` and `collapse` steps. Check their [documentation](https://flair.readthedocs.io/en/stable/) for more information.
    - Additionally FLAIR recommends to split bed files that are larger than 1GB into smaller chromosome specific fractions and concatenate them later again. You can configure the breakpoint if the runtime on your system requires smaller or allows larger files.
  - IsoQuant needs to know the type of the reads. `--complete_genedb` is recommended when using official reference annotations. Check their [documentation](https://ablab.github.io/IsoQuant/) for more information.
  - IsoTools is configured independently for each run. Each version can pass a dict of settings to `add_sample_from_bam`, `add_qc_metrics`, `write_gtf` and `training_sequences`. The only required parameter is `query` for `write_gtf`. Check their [documentation](https://isotools.readthedocs.io/en/latest/) for more information.
  - StringTie requires the path to its executable. More parameters for both the `run` and the `merge` step can be passed. Check their [documentation](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) for more information.
  - STAR only runs when short reads are available. You can configure extra parameters for the `index` and the `align` step. Check their (extensive) [documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more information.
  - SQANTI3 requires multiple things
    - the path to the cloned SQANTI3 repo to access their scripts
    - the path to a polyA motif list. You can use the one provided in the repo.
    - additional parameters for the `qc` step. Check their [documentation](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control) for more information.
    - multiple heatmaps are created based on the sqanti output of all tools and sample groups, which are defined with `plot_groups`
      - Specify which tools and groups you want to have for each plot group and optionally also specify overwrite names for both.
      - Specify if you want to have titles in the plots. Useful for presentations, but usually not needed for publications.
      - set the desired DPI
      - specify the color maps for TSS and PAS plots
- If snakemake has trouble to find conda, uncomment `path_to_conda` and specifiy the path to the conda executable. This might or might not help.

# Output Files

*All paths are relative to the `snakemake-pipeline/` folder.*

## Logs

Each rule produces at least one log file, which can be found in `logs/`. Errors from some of the python scripts are hard to capture, so they might land either in stdout/stderr or in log files from your cluster environment depending on your snakemake configuration.

## Results

Each tool produces one gtf file for each group of samples. These can be found in `results/<tool>/transcriptome/<group>.gtf`. A gzipped file with a tabix index is also available.

If you're interested in some of the tool specific output files, you can find them in the respective tool folder in `results/<tool>/`.

The SQANTI3 results of each tool are stored in `results/sqanti/<tool>/<group>`. You'll find the SQANTI3 report in HTML and PDF format, as well as the `<group>_classification.txt` file.

## Plots

All plots are in `results/plots`.

The overlap between tools is in `results/plots/upset/`. `all/` contains all combinations, while `filtered/` excludes the combinations with only one tool. Both folders also have a `mean.png` file, that uses the average overlap over all groups.

The heatmaps are in `results/plots/sqanti/`. One folder for each plot group, which you defined in the config. For each `CAGE_support`, `TSS_ratio`, `polyA_motif` and `polyA_site` at least 20 different heatmaps are created. Additionally three `transcript_counts` bar plot show the distribution of transcripts in the SQANTI3 categories and subcategories. And last there is `stats.tsv`, which outputs the correlation between between CAGE support and TSS ratio, as well as the correlation between polyA motif and polyA site.
