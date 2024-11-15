#!/usr/bin/env python3

import pandas as pd
import re

def transform_value(value):
    # Split the value into components
    if value != "NA|NA|NA":
        parts = value.split('|')
        if len(parts) == 3:
            # Reorder and format the components into the desired format
            nuc_identity = parts[0][:-2]  # Extract the number from '98NT'
            aa_identity = parts[1][:-2]   # Extract the number from '98AA'
            coverage = parts[2]           # Extract the coverage number
            # Format the new string
            return f'[{coverage}NT/{aa_identity}AA/{nuc_identity}]G'
    else:
        return ""
    return value

def clean_and_format_centar_dfs(centar_df):
    '''If Centar was run get info to add to the dataframe.'''
    cols_to_transform = [x for x in centar_df.columns if '[%Nuc_Identity' in x ]
    for col in cols_to_transform:
        centar_df[col] = centar_df[col].apply(transform_value)
    #drop presence/absence columns
    columns_to_drop = [col for col in centar_df.columns if 'presence'  in col ]
    clean_centar_df = centar_df.drop(columns=columns_to_drop)
    # Remove the substring from all column headers
    clean_centar_df.rename(columns=lambda x: re.sub(r'\[%Nuc_Identity \| %AA_Identity \| %Coverage\]', '', x).strip(), inplace=True)
    clean_centar_df.rename(columns=lambda x: re.sub(r'Diffbase_', '', x).strip(), inplace=True)
    #separate dataframes
    RB_type = [ "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info" ]
    RB_type_col = [col for col in clean_centar_df.columns if any(substring in col for substring in RB_type) ]
    RB_type_len = len(RB_type_col)
    A_B_Tox = [ "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB"]
    A_B_Tox_col = [col for col in clean_centar_df.columns if any(substring in col for substring in A_B_Tox) ]
    A_B_Tox_len = len(A_B_Tox_col)
    other_Tox = [ "tcdC", "tcdR", "tcdE", "cdtA", "cdtB", "cdtR", "cdtAB1", "cdtAB2", "non-tox" ]
    other_Tox_col = [col for col in clean_centar_df.columns if any(substring in col for substring in other_Tox) ]
    other_Tox_len = len(other_Tox_col)
    mutants = [ 'gyr','dac','feo','fur','gdp','gly','hem','hsm','Isc','mur', 'mur','nifJ','PNim','rpo','sda','thi','Van','mutations' ]
    mutations_col = [col for col in clean_centar_df.columns if any(substring in col for substring in mutants) ]
    mutant_len = len(mutations_col)
    existing_columns_in_order = "MLST Clade" + A_B_Tox_col + other_Tox_col + mutations_col + RB_type_col
    ordered_centar_df = clean_centar_df[existing_columns_in_order]
    #df_cleaned = df.replace('', np.nan).dropna(axis=1, how='all')
    return ordered_centar_df, A_B_Tox_len, other_Tox_len, RB_type_len

def create_centar_combined_df(directory, sample_name):
    '''If Centar was run get info to add to the dataframe.'''
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    # create file names
    centar_summary = directory + "/" + sample_name + "/CENTAR/" + sample_name + "_centar_output.tsv"
    #clean up the dataframe
    centar_df = pd.read_csv(centar_summary, sep='\t', header=0)
    return centar_df


######################################## ShigaPass functions ##############################################
def create_shiga_df(directory, sample_name, shiga_df):
    '''If Shigapass was run get info to add to the dataframe.'''
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    # create file names
    shiga_summary = directory + "/ANI/" + sample_name + "_ShigaPass_summary.csv"
    # Create a dictionary to store row information
    row_data = { "WGS_ID": sample_name, "ShigaPass_Organism": ""}
    row_data["WGS_ID"] = sample_name
    with open(shiga_summary) as shiga_file:
        for line in shiga_file.readlines()[1:]:  # Skip the first line
            if line.split(";")[9] == 'Not Shigella/EIEC\n':
                row_data["ShigaPass_Organism"] = ""
            else:
                row_data["ShigaPass_Organism"] = line.split(";")[7]
            # Convert the row data into a DataFrame and concatenate with the main DataFrame
            shiga_df = pd.concat([shiga_df, pd.DataFrame([row_data])], ignore_index=True)
    # Define the mapping of short strings to longer strings
    mapping_dict = {'SB': 'Shigella boydii',
                    'SD': 'Shigella dysenteriae',
                    'SS': 'Shigella sonnei',
                    'SF1-5': 'Shigella flexneri'}
    # Apply the mapping using map()
    shiga_df['ShigaPass_Organism'] = shiga_df['ShigaPass_Organism'].replace(mapping_dict)
    return shiga_df

def double_check_taxa_id(shiga_df, phx_df):
    # Merge the DataFrames on 'WGS_ID'
    merged_df = pd.merge(phx_df, shiga_df, on='WGS_ID', how='left')
    # Identify the position of the insertion point
    insert_position = merged_df.columns.get_loc("FastANI_Organism")
    # Reorder the columns: place the new columns at the desired position
    columns = list(merged_df.columns)
    # Reorder columns to insert the new columns between 'Column_A' and 'Column_B'
    new_columns = ['ShigaPass_Organism']
    # Reorder columns: place the new columns between 'Column_A' and 'Column_B'
    columns_reordered = (
        columns[:insert_position] +  # Columns before the insertion point
        new_columns +                # New columns to be inserted
        [col for col in columns if col not in new_columns and col not in columns[:insert_position]]) # Remaining columns
    # Reorder the merged DataFrame columns
    merged_df = merged_df[columns_reordered]
    # Apply the custom function to fill the Taxa_ID column
    merged_df['Final_Taxa_ID'] = merged_df.apply(fill_taxa_id, axis=1)
    return merged_df

#def main():
#    directory = "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Jill_DIR/PHX_v2/v2.2.0-dev/centar/cdc_centar_newer"
#    sample_names = [ "2022GL-00907", "2022GL-00947", "2022GL-01162" ]
#    centar_dfs = []
#    for sample_name in sample_names:
#        centar_df = create_centar_combined_df(directory, sample_name)
#        centar_dfs.append(centar_df)
#    full_centar_df = pd.concat(centar_dfs, ignore_index=True)
#    ordered_centar_df = clean_and_format_centar_dfs(full_centar_df)
