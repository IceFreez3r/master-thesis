import pandas as pd
import numpy as np
from upsetplot import from_memberships, plot
import matplotlib.pyplot as plt

import logging

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])

logger = logging.getLogger("overlap")

# Import all dataframes from tsv
dfs = [pd.read_csv(file, sep="\t") for file in snakemake.input]
# Concat
df = pd.concat(dfs).reset_index(drop=True)
# Average
df = df.groupby("tools").mean().reset_index()
# Parse tools from string
df["tools"] = df["tools"].apply(
    lambda s: [e.replace('(', '').replace(')', '').replace(' ', '').replace("'", '')
               for e in s.split(",")]).apply(lambda x: [e for e in x if len(e) > 0])
# Sort
vlen = np.vectorize(len)
df = df.sort_values("tools", key=lambda x: vlen(x))
print(df["tools"].to_list())

logger.info(f"Creating upset plot for all tissues")
upset = from_memberships(df["tools"].to_list(), data=df["unique"].to_list())
plot(upset, sort_categories_by='input', orientation='horizontal', element_size=30)
plt.savefig(snakemake.output.upset)
plt.close()

# Exclude single tool
df_filtered = df[df["tools"].apply(len) > 1]
upset = from_memberships(df_filtered["tools"].to_list(), data=df_filtered["unique"].to_list())
plot(upset, sort_categories_by='input', orientation='horizontal', element_size=30)
plt.savefig(snakemake.output.upset_filtered)
plt.close()

logger.info("Done")
