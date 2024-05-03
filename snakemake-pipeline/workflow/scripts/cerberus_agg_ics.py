import subprocess
import pandas as pd
import os


df = pd.DataFrame({'BED file path': [snakemake.input.ref], 'Reference': [True], 'Source name': ['ENCODE']})
for ics in snakemake.input.ics:
    df = pd.concat([df, pd.DataFrame({'BED file path': [ics], 'Reference': [False], 'Source name': ['TALON']})])

with open(snakemake.log[0], 'w', buffering=1) as f:
    if snakemake.params.limit is None:
        df.to_csv(snakemake.output.config, index=False, header=False)
    else:
        # Split the dataframe into chunks of size limit
        chunks = [df[i:i+snakemake.params.limit] for i in range(0, len(df), snakemake.params.limit)]
        f.write(f'Aggregating in {len(chunks)} chunks of length {snakemake.params.limit}\n')
        for i, chunk in enumerate(chunks):
            # prepend the output of each chunk to the next one
            if i != 0:
                # Source "cerberus" to indicate that
                # > it shouldn't modify any of the already-aggregated sources that are already in the ICs file
                chunk = pd.concat([pd.DataFrame({'BED file path': [f'{snakemake.output.tsv}.{i-1}'], 'Reference': [False], 'Source name': ['cerberus']}), chunk])
            if i == len(chunks) - 1:
                chunk.to_csv(snakemake.output.config, index=False, header=False)
            else:
                config = f'{snakemake.output.config}.{i}'
                output = f'{snakemake.output.tsv}.{i}'
                if os.path.exists(output):
                    f.write(f'{output} already exists, skipping\n')
                    continue
                chunk.to_csv(config, index=False, header=False)
                f.write(f'cerberus agg_ics --input {config} -o {output}\n')
                subprocess.run(f'cerberus agg_ics --input {config} -o {output}', shell=True, stdout=f, stderr=f, check=True)

    f.write(f'cerberus agg_ics --input {snakemake.output.config} -o {snakemake.output.tsv}\n')
    subprocess.run(f'cerberus agg_ics --input {snakemake.output.config} -o {snakemake.output.tsv}', shell=True, stdout=f, stderr=f, check=True)
