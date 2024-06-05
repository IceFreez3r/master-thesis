import requests
import os
import pandas as pd

sample_df = pd.read_csv("/project/hfa_work/ENCODE/data/reads/metadata_tissue.tsv", sep="\t")

encode_url = "https://www.encodeproject.org"

data_dir = "/project/hfa_work/ENCODE/data/talon_gtf/"

gtf_metadata = pd.DataFrame(columns=['sample ID', 'gtf'])

for i, sample_row in sample_df.iterrows():
    sample_id = sample_row['sample ID']
    print(f"Processing sample {sample_id}")
    # get JSON for each sample (sample ID) from ENCODE
    url = f"{encode_url}/experiments/{sample_id}/?format=json"
    response = requests.get(url)
    response_json = response.json()
    print("Got JSON response from ENCODE.")
    # check files.file_type for 'gtf'
    for file in response_json['files']:
        if file['file_type'] == 'gtf':
            gtf_url = file['href']
            break
    gtf_url = encode_url + gtf_url
    filename = gtf_url.split("/")[-1]
    if not os.path.exists(f"{data_dir}{filename}"):
        print(f"Downloading GTF file for sample {sample_id}")
        os.system(f"wget -P {data_dir} {gtf_url}")
        print(f"Downloaded GTF file for sample {sample_id}")
    else:
        print(f"GTF file for sample {sample_id} already exists.")
    gtf_metadata = pd.concat([gtf_metadata, pd.DataFrame({'sample ID': sample_id, 'gtf': f"{data_dir}{filename}"}, index=[0])], ignore_index=True)

gtf_metadata.to_csv(f"{data_dir}gtf_metadata.tsv", sep="\t", index=False)
