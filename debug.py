class Foo(object):
    pass


import pandas as pd
#samples_table = pd.read_csv("samples.csv").set_index("sample", drop=False)
units = pd.read_table('config/samples.tsv').set_index("sample_name",drop=False)

# make a variable, wildcards, which is an object of that class

wildcards = Foo()

# now, if you want to test specific values for different wildcards you can do like:
wildcards.sample = "s002"

# print that value:
wildcards.sample


def fq_dict_from_sample(wildcards):
  return {
    "fq1": samples_table.loc[wildcards.sample, "fastq1"], 
    "fq2": samples_table.loc[wildcards.sample, "fastq2"]
  }

# make sure our spoofed wildcards variable is set:
wildcards.sample = "s003"

# Now see what that function returns on that wildcards input:
fq_dict_from_sample(wildcards)

