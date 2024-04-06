import os
import pandas as pd
import config

raw_data_dir = f"{config.base_output_dir}/Raw_data/Illumina-ShortReads"
data_files = os.listdir(raw_data_dir)
pd.DataFrame(list(set(["_".join(d.split("_")[:-3]) for d in data_files]))).to_csv("acc_list.csv",index=None,header=None)