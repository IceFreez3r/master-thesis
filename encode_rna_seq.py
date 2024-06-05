import requests
import os
import pandas as pd

sample_df = pd.read_csv("/project/hfa_work/ENCODE/data/reads/metadata_tissue.tsv", sep="\t")

encode_url = "https://www.encodeproject.org"

data_dir = "/project/hfa_work/ENCODE/data/rna_seq/"

rnaseq_metadata = pd.DataFrame(columns=['sample ID', 'BAM'])

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
    for file in response_json['files']:
        if file['file_type'] == 'bam' and file['output_type'] == 'alignments':
            BAM_url = file['href']
            break
    else:
        print(f"No BAM file found for sample {sample_id}")
        continue
    BAM_url = encode_url + BAM_url
    filename = BAM_url.split("/")[-1]
    if not os.path.exists(f"{data_dir}{filename}"):
        print(f"Downloading BAM file for sample {sample_id}")
        os.system(f"wget -P {data_dir} {BAM_url}")
        print(f"Downloaded BAM file for sample {sample_id}")
    else:
        print(f"BAM file for sample {sample_id} already exists.")
    rnaseq_metadata = pd.concat([rnaseq_metadata, pd.DataFrame({'sample ID': sample_id, 'BAM': f"{data_dir}{filename}"}, index=[0])], ignore_index=True)

rnaseq_metadata.to_csv(f"{data_dir}rnaseq_metadata.tsv", sep="\t", index=False)
