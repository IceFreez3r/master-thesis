import pandas as pd
import glob
import os

def get_sample_df(config):
    sample_df = pd.read_csv(config['sample_table'], sep="\t")
    return sample_df

def get_samples(config):
    sample_df = get_sample_df(config)
    samples = sample_df["sample ID"].tolist()
    return samples

def get_tissues(config):
    sample_df = get_sample_df(config)
    tissues = sample_df["group"].unique().tolist()
    CAGE_dir = config["CAGE_dir"]
    tissues = [tissue for tissue in tissues if glob.glob(f"{CAGE_dir}/*{tissue}*.bed.gz")]
    return tissues

def get_experiments(config):
    sample_df = get_sample_df(config)
    experiments = sample_df["file"].tolist()
    experiments = [experiment.split("/")[-1].split(".")[0] for experiment in experiments]
    return experiments


def get_experiment_for_sample(config, sample):
    sample_df = get_sample_df(config)
    experiment = sample_df[sample_df["sample ID"] == sample]["file"].values[0]
    experiment = experiment.split("/")[-1].split(".")[0]
    return experiment

def get_tissue_for_sample(config, sample):
    sample_df = get_sample_df(config)
    tissue = sample_df[sample_df["sample ID"] == sample]["group"].values[0]
    return tissue

def get_alignment_for_sample(config, sample):
    return os.path.join(config['alignment_dir'], get_experiment_for_sample(config, sample) + '_aligned.bam')


def get_samples_for_tissue(config, tissue):
    sample_df = get_sample_df(config)
    samples = sample_df[sample_df["group"] == tissue]["sample ID"].tolist()
    return samples

def get_experiments_for_tissue(config, tissue):
    sample_df = get_sample_df(config)
    experiments = sample_df[sample_df["group"] == tissue]["file"].tolist()
    experiments = [experiment.split("/")[-1].split(".")[0] for experiment in experiments]
    return experiments
