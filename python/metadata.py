#%%
import pandas as pd
import re
#%%
metadata = pd.read_excel('data/EI_key_rainbow.xlsx')
#%%
def extract_ipas_part(filename):
    pattern = r'(IPAS.*?)\.htrms'
    match = re.search(pattern, filename)
    return match.group(1) if match else None

# Apply the function to the 'FileName' column
metadata['PLateNumber'] = metadata['PLateNumber'].apply(extract_ipas_part)

# create group variable based on Is_Case, if Is_Case is 1, group is 'Case', and if 0 as  'Control', other as 'Unknown'
metadata['group'] = metadata['Is_Case'].map({1: 'Case', 0: 'Control'})

metadata
#%%
metadata['Cancer Type'].value_counts()
#%%
metadata.to_csv('data/metadata.csv', index=False)
#%%

#%%
