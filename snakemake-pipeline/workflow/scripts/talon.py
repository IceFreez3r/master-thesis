import os
import pandas as pd

use_tmpdir = snakemake.params.get("use_tmpdir", False)

samples = snakemake.params.samples
sample_df = pd.read_csv(snakemake.config["sample_table"], sep="\t")


def talon_config():
    config = []
    for i, sample in enumerate(samples):
        config.append(
            {
                "dataset_name": sample,
                "sample_description": sample_df.iloc[i]["Description"],
                "platform": snakemake.params.sequencing_platform,
                "sam": snakemake.input.sample_sams[i],
            }
        )
        if use_tmpdir:
            # copy sam file to tmpdirectory and change the path in the config
            new_sam = os.path.join(snakemake.resources.tmpdir, f"{sample}.sam")
            os.system(f"cp {snakemake.input.sample_sams[i]} {new_sam}")
            config[-1]["sam"] = new_sam
    config = pd.DataFrame(config)
    config.to_csv(snakemake.output.config, index=False, header=False)


print("Creating TALON config file...")

talon_config()

print("Config file created. Starting TALON run...")

os.system(
    f"talon --db {snakemake.input.db} --build {snakemake.params.genome_name} --threads {snakemake.threads} --o {snakemake.params.output_prefix} --f {snakemake.output.config} --tmpDir {snakemake.resources.tmpdir}/talon_tmp > {snakemake.log} 2>&1"
)

print("TALON run complete")
