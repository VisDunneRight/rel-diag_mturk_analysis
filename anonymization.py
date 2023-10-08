# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Anonymization of user data
#
# Script that takes as in put the raw data from AMT and anonymizes the data:
# Replaces the worker ID and removes 'assignment_id', 'hit_id', and 'qualification_score'.
# This script is for internal use and will not be shared.
#

# %% [markdown] jupyter={"outputs_hidden": false}
# ## Load packages

# %%
import pandas as pd
import numpy as np
from IPython.display import display

# Set Jupyter and Pandas to show 3 decimal places, does not work for lists of numbers
# %precision 3
pd.options.display.float_format = '{:,.3f}'.format
np.set_printoptions(precision=3)

# %% [markdown]
# ## File names

# %%
originalfilename = 'data/users-table.csv'                                  # original data: will not be shared
newfilenameanonymized = 'data/users-table-anonymized.csv'                  # anonymized data: will be used later
internalencodingfilename = 'data/users-table-original-worker_ids.csv'      # for internal use: store the transformation from anonymized to original worker_ids

# %% [markdown]
# ## Creating two new files
#
# Loading the original data, anonymizing it, saving it and also the worker-id encoding
#

# %%
df = pd.read_csv(originalfilename)
print(f'Original data: {originalfilename}')
display(df)

# --- Anonymize the workers with randomized categories (https://pandas.pydata.org/docs/user_guide/categorical.html#working-with-categories)
c = df.worker_id.astype('category')
d = dict(enumerate(c.cat.categories))

import random
from random import shuffle
random.seed(42)# so we get the same result each time

keys = list(d.keys())
shuffle(keys)
d = dict(zip(d.values(), keys))             # randomize the category assignment s.t. categorical numbers are identical with alphabetical order

with open(internalencodingfilename, 'w') as f:      # save the category encoding
    for key in d.keys():
        f.write("%s,%s\n"%(key,d[key]))

df["worker_id"] = df["worker_id"].map(d)    # replace worker_ids with randomized categories
# df.drop(['assignment_id', 'hit_id', 'qualification_score'], axis=1, inplace=True)     # keep original schema
df.loc[df['assignment_id'] != '', 'assignment_id'] = ''
df.loc[df['hit_id'] != '', 'hit_id'] = ''
df.loc[df['tutorial_time'] != '', 'tutorial_time'] = ''# Blank tutorial times because they were collected erroneously

# --- store anonymized data and print
df.to_csv(newfilenameanonymized, index=False)
print(f'New data: {newfilenameanonymized}')
display(df)

# %% [markdown] jupyter={"outputs_hidden": false}
# # Evaluate random success rate

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
# What is the probability that a user randomly guessing answers k/n questions correctly? Binomial distribution
from scipy.stats import binom

# 50% correct 16/32
k=16
n=32
p=1/4
print(1-binom.cdf(k-1, n, p))

# 66.6% correct 24/32
k=24
n=32
p=1/4
print(1-binom.cdf(k-1, n, p))

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
# What is the probability that among n participants, there are two who have the same treatment?
n = 50
p = 1 - np.exp(-n*n/(2*2520*2520))
print(p*5000)
