import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
import itertools
import duckdb
from upsetplot import from_memberships, plot
import logging

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])

logger = logging.getLogger("isotools")
# Disable future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

tools = snakemake.params.tools
tss_error = snakemake.params.tss_error
pas_error = snakemake.params.pas_error
junction_error = snakemake.params.junction_error

def get_tabix(tool):
    return pysam.TabixFile(snakemake.input[tool])

def tool_transcriptome(tool):
    transcripts = []
    exon_list = []
    exon_count = {}
    tabix = get_tabix(tool)
    for row in tabix.fetch(parser=pysam.asGTF()):
        if row.feature == "transcript":
            transcripts.append((row.transcript_id, row.contig, int(row.start), int(row.end), row.strand))
        elif row.feature == "exon":
            exon_count.setdefault(row.transcript_id, 0)
            exon_list.extend([
                (row.transcript_id, 2 * exon_count[row.transcript_id], int(row.start)),
                (row.transcript_id, 2 * exon_count[row.transcript_id] + 1, int(row.end))])
            exon_count[row.transcript_id] += 1

    transcriptome = pd.DataFrame(transcripts, columns=["transcript_id", "chr", "start", "end", "strand"])
    transcriptome["min_start"] = np.where(transcriptome["strand"] == '+', transcriptome["start"] - tss_error, transcriptome["start"] - pas_error)
    transcriptome["max_start"] = np.where(transcriptome["strand"] == '+', transcriptome["start"] + tss_error, transcriptome["start"] + pas_error)
    transcriptome["min_end"] = np.where(transcriptome["strand"] == '+', transcriptome["end"] - pas_error, transcriptome["end"] - tss_error)
    transcriptome["max_end"] = np.where(transcriptome["strand"] == '+', transcriptome["end"] + pas_error, transcriptome["end"] + tss_error)

    exons = pd.DataFrame(exon_list, columns=["transcript_id", "index", "position"])
    transcriptome.set_index("transcript_id", inplace=True)
    transcriptome["max_splice_site_index"] = exons.groupby("transcript_id").max()["index"]
    transcriptome.reset_index(inplace=True, names=["transcript_id"])

    return transcriptome, exons


def add_transcript_overlap_duckdb(tool, tool_overlap, overlap_exons, remove_duplicates=True):
    transcriptome, new_exons = tool_transcriptome(tool)

    if tool_overlap.empty:
        tool_overlap = pd.concat([tool_overlap, transcriptome])
        tool_overlap[tool] = transcriptome["transcript_id"]
        tool_overlap["id"] = tool + transcriptome["transcript_id"]
        tool_overlap.drop(columns=["transcript_id"], inplace=True)

        new_exons['id'] = tool + new_exons["transcript_id"]
        new_exons.drop(columns=["transcript_id"], inplace=True)
        overlap_exons = pd.concat([overlap_exons, new_exons])

        return tool_overlap, overlap_exons

    query = f"""
SELECT
    tool_ol.id,
    t.transcript_id
FROM tool_overlap AS tool_ol
JOIN transcriptome t ON
    tool_ol.chr = t.chr AND
    tool_ol.start >= t.min_start AND
    tool_ol.start <= t.max_start AND
    tool_ol.end >= t.min_end AND
    tool_ol.end <= t.max_end AND
    tool_ol.strand = t.strand AND
    tool_ol.max_splice_site_index = t.max_splice_site_index
JOIN new_exons ne ON
    ne.transcript_id = t.transcript_id
JOIN overlap_exons oe ON
    oe.id = tool_ol.id AND
    oe.index = ne.index AND
    (
        ABS(oe.position - ne.position) <= {junction_error} OR
        oe.index = 0 OR
        oe.index = t.max_splice_site_index
    )
GROUP BY
    tool_ol.id,
    t.transcript_id
HAVING
    COUNT(*) = MAX(t.max_splice_site_index) + 1
               """

    matches = duckdb.sql(query).df()

    # Happens when multiple transcripts from the new tool match the same transcript from the previous tools
    # or when a transcript from the new tool matches multiple transcripts from the previous tools
    if remove_duplicates:
        logger.info("with duplicates:", matches.shape)
        matches = matches.drop_duplicates(subset=["id"]).drop_duplicates(subset=["transcript_id"])
        logger.info("without duplicates", matches.shape)

    # Add matches to the df
    tool_overlap = tool_overlap.merge(matches, on="id", how="left")
    tool_overlap[tool] = tool_overlap["transcript_id"]
    tool_overlap.drop(columns=["transcript_id"], inplace=True)

    # Add mismatches to the df
    transcriptome_filtered = transcriptome.set_index("transcript_id").drop(matches["transcript_id"].to_list()).reset_index(names=["transcript_id"])
    transcriptome_filtered[tool] = transcriptome_filtered["transcript_id"]
    transcriptome_filtered["id"] = tool + transcriptome_filtered["transcript_id"]
    transcriptome_filtered.drop(columns=["transcript_id"], inplace=True)
    tool_overlap = pd.concat([tool_overlap, transcriptome_filtered])

    # Add exons for mismatches
    new_exons['id'] = tool + new_exons["transcript_id"]
    new_exons.drop(columns=["transcript_id"], inplace=True)
    # Filter new_exons to only contain exons of transcriptome_filtered
    new_exons = new_exons[new_exons["id"].isin(transcriptome_filtered["id"])]
    overlap_exons = pd.concat([overlap_exons, new_exons])

    return tool_overlap, overlap_exons

overlap = pd.DataFrame(columns=["chr", "start", "end", "strand", "max_splice_site_index", *tools])
overlap_exons = pd.DataFrame(columns=["id", "index", "position"])

for tool in tools:
    logger.info("Processing", tool)
    overlap, overlap_exons = add_transcript_overlap_duckdb(tool, overlap, overlap_exons)

# iterate over all combination of tools and count the number of transcripts that are shared
combinations = itertools.chain.from_iterable(itertools.combinations(tools, r) for r in range(1, len(tools) + 1))

counts = []
for combination in combinations:
    inverted_combination = [tool for tool in tools if tool not in combination]
    counts.append((combination,
                   overlap[[*combination]].dropna().shape[0],
                   overlap[(overlap[inverted_combination].isna().all(axis=1)) & (overlap[[*combination]].notna().all(axis=1))].shape[0]))

counts = pd.DataFrame(counts, columns=["tools", "shared", "unique"])

upset = from_memberships(counts["tools"].to_list(), data=counts["unique"].to_list())
plot(upset, orientation='horizontal', element_size=30)
plt.savefig(snakemake.output.upset)
plt.show()

# remove stringtie and flair unique entries
counts_filtered = counts[(counts["tools"] != ("flair",)) & (counts["tools"] != ("stringtie",))]
upset = from_memberships(counts_filtered["tools"].to_list(), data=counts_filtered["unique"].to_list())
plot(upset, orientation='horizontal', element_size=30)
plt.savefig(snakemake.output.upset_filtered)
plt.show()
