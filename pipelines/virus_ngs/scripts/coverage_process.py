
from pathlib import Path
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# plt.style.use('ggplot')

### constants
win_size = 10
cov_cutoffs = [50, 100]

config_df = pd.read_csv(snakemake.input.sample_table)

# coverage data
with open(snakemake.output.where_at, "w") as f:
    full_cov_df = None
    for coverage_file in snakemake.input.coverage_files:
        coverage_file = Path(coverage_file)
        sample = coverage_file.name.split(".")[0]
        header = ['chrom', 'pos', 'coverage']
        f.write("processing: {}\n".format(sample))
        
        df = pd.read_csv(coverage_file, sep="\t", names=header)
        df['sample'] = sample
        # df['bioproject'] = df['sample'].map(mapping_dict)
        
        # mask first and last 100bp
        df = df.head(-100).tail(-100)

        # rolling window mean
        df['win_cov'] = df['coverage'].rolling(window=win_size, min_periods=2, center=True).mean()

        # clean
        # columns = ["sample", "bioproject", "pos", "coverage", "win_cov"]
        columns = ["sample", "pos", "coverage", "win_cov"]
        df = df[columns]

        # can do plotting here...

        # initialize the coverage results df 
        # cov_df = pd.DataFrame(data={"sample":[sample], "bioproject":[mapping_dict[sample]]})
        cov_df = pd.DataFrame(data={"sample":[sample]})

        # mean, median coverage
        cov_df['cov_mean'] = df['coverage'].mean().round(decimals=2)
        cov_df['cov_median'] = df['coverage'].median()

        for cov_cutoff in cov_cutoffs:
            column_cov = "genome_cov_{}X".format(cov_cutoff)
            cov_df[column_cov] = len(df[df['coverage'] >= cov_cutoff]) / df.shape[0]

            
        if full_cov_df is None:
            full_cov_df = cov_df.copy()
        else:
            full_cov_df = pd.concat([full_cov_df, cov_df], ignore_index=True)

# merge with config_df
config_df = pd.merge(config_df, full_cov_df, on='sample')

# write coverage summary to file
# full_cov_df.to_csv(snakemake.output.coverage_QC_file, index=False)
config_df.to_csv(snakemake.output.coverage_QC_file, index=False)
