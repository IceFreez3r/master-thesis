import os
import pandas as pd

use_tmpdir = snakemake.params.get("use_tmpdir", False)

sample_df = pd.read_csv(snakemake.config["sample_table"], sep="\t")
samples = sample_df["file"].apply(lambda x: x.split("/")[-1].split(".")[0])

def talon_config():
    config = []
    for sample in samples:
        sample_index = samples[samples == sample].index[0]
        config.append({
            'dataset_name': sample,
            'sample_description': sample_df.iloc[sample_index]['Description'],
            'platform': snakemake.params['sequencing_platform'],
            'sam': snakemake.input['sample_sams'][sample_index]
        })
        if use_tmpdir:
            # copy sam file to tmpdirectory and change the path in the config
            new_sam = os.path.join(snakemake.tmpdir, f"{sample}.sam")
            os.system(f"cp {snakemake.input['sample_sams'][sample_index]} {new_sam}")
            config[-1]['sam'] = new_sam
    config = pd.DataFrame(config)
    config.to_csv(snakemake.output['config'], index=False, header=False)

talon_config()

os.system(f"talon --db {snakemake.input['db']} --build {snakemake.params['genome_name']} --threads {snakemake.threads} --o {snakemake.params['output_prefix']} --f {snakemake.output['config']} > {snakemake.log} 2>&1")

print("TALON run complete")
