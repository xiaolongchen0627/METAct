import glob

import pandas as pd
from snakemake.utils import validate

samples = ( 
pd.read_table("samples.tsv")
.set_index("sample",drop=False))