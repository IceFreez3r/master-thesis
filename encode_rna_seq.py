import requests
import os
import pandas as pd

sample_df = pd.read_csv("/path/to/data/reads/metadata_tissue.tsv", sep="\t")

encode_url = "https://www.encodeproject.org"

data_dir = "/path/to/data/rna_seq/"

rnaseq_metadata_bam = pd.DataFrame(columns=['sample ID', 'file'])
rnaseq_metadata_fastq = pd.DataFrame(columns=['sample ID', 'file'])

for i, sample_row in sample_df.iterrows():
    sample_id = sample_row['sample ID']
    print(f"Processing sample {sample_id}")
    # get JSON for each sample (sample ID) from ENCODE
    url = f"{encode_url}/experiments/{sample_id}/?format=json"
    response = requests.get(url)
    response_json = response.json()
    print("Got JSON response from ENCODE.")

    # get biosample id
    biosample = response_json['replicates'][0]['library']['biosample']['accession']
    print(f"Biosample: {biosample}")
    url = f"{encode_url}/search/?type=Experiment&searchTerm={biosample}&format=json"
    response = requests.get(url)
    response_json = response.json()
    print("Got JSON response from ENCODE.")
    for experiment in response_json['@graph']:
        if experiment['assay_term_name'] == 'RNA-seq':
            rna_seq_id = experiment['accession']
            break
    else:
        print(f"No RNA-seq experiment found for biosample {biosample}")
        continue
    print(f"RNA-seq ID: {rna_seq_id}")

    # get rna-seq reads
    url = f"{encode_url}/experiments/{rna_seq_id}/?format=json"
    response = requests.get(url)
    response_json = response.json()
    print("Got JSON response from ENCODE.")
    if (response_json['biosample_summary'] != sample_row['Biosample summary']):
        print(f"Biosample summary does not match for sample {sample_id}:\n{response_json['biosample_summary']}\nvs\n{sample_row['Biosample summary']}")
        continue

    BAM_url, fastq_url = None, None
    for file in response_json['files']:
        if file['file_type'] == 'bam' and file['output_type'] == 'alignments':
            BAM_url = file['href']
        if file['file_type'] == 'fastq' and file['output_type'] == 'reads':
            fastq_url = file['href']
    if BAM_url is None:
        print(f"No BAM file found for sample {sample_id}")
        continue
    if fastq_url is None:
        print(f"No fastq file found for sample {sample_id}")
        continue

    BAM_url = encode_url + BAM_url
    filename = BAM_url.split("/")[-1]
    if not os.path.exists(os.path.join(data_dir, filename)):
        print(f"Downloading BAM file for sample {sample_id}")
        os.system(f"wget -P {data_dir} {BAM_url}")
        print(f"Downloaded BAM file for sample {sample_id}")
    else:
        print(f"BAM file for sample {sample_id} already exists.")
    rnaseq_metadata_bam = pd.concat([rnaseq_metadata_bam, pd.DataFrame({'sample ID': sample_id, 'file': f"{data_dir}{filename}"}, index=[0])], ignore_index=True)

    fastq_url = encode_url + fastq_url
    filename = fastq_url.split("/")[-1]
    if not os.path.exists(os.path.join(data_dir, filename)):
        print(f"Downloading fastq file for sample {sample_id}")
        os.system(f"wget -P {data_dir} {fastq_url}")
        print(f"Downloaded fastq file for sample {sample_id}")
    else:
        print(f"Fastq file for sample {sample_id} already exists.")
    rnaseq_metadata_fastq = pd.concat([rnaseq_metadata_fastq, pd.DataFrame({'sample ID': sample_id, 'file': f"{data_dir}{filename}"}, index=[0])], ignore_index=True)

rnaseq_metadata_bam.to_csv(os.path.join(data_dir, "rnaseq_metadata_bam.tsv"), sep="\t", index=False)
rnaseq_metadata_fastq.to_csv(os.path.join(data_dir, "rnaseq_metadata_fastq.tsv"), sep="\t", index=False)
