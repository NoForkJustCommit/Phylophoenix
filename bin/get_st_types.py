#!/usr/bin/env python3

import pandas as pd
import argparse
import csv

##Makes a new samplesheet for assemblies based on their st type. Expects that assemblies are in an output PHoeNIx folder
##Usage: >python get_st_types.py -s GRiPHin_samplesheet.csv -g GRiPHin_Report.xlsx
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--samplesheet', default=None, required=False, dest='samplesheet', help='GRiPHin samplesheet of sample,directory in csv format.')
    parser.add_argument('-g', '--griphin_report', required=False, dest='griphin_report', help='A griphin excel report.')
    return parser.parse_args()

def get_st_groups(griphin_report):
    """Get unique ST groups from griphin report."""
    df = pd.read_excel(griphin_report, header=1)
    df.dropna(subset=['Primary_MLST'], inplace=True) #drop rows where there is na in the MLST column
    clean_df = df[~df.Primary_MLST.str.contains("-")] #remove samples that do not have MLST
    novel_df = clean_df[clean_df["Primary_MLST"].str.contains("Novel_allele")] # if there is novel allele in the primary mlst column separate them out
    indexes_to_drop = list(novel_df.index) # get indexes that have "Novel_allele" in "Primary_MLST" column
    clean_df = clean_df.drop(index=indexes_to_drop) # drop the rows that have "Novel_allele" in "Primary_MLST" column
    list_of_sts = clean_df["Primary_MLST"].unique() # get unique mlsts from what is left
    #get sample names per st type into dictionary
    st_dict = {} #create an empty dictionary
    for seq_type in list_of_sts:
        sample_list = list(clean_df.loc[clean_df["Primary_MLST"] == seq_type, 'WGS_ID']) # get the sample names per ST type
        if len(sample_list) > 1: # Do not include cases where the there is only one sample for a given ST type
            st_dict[seq_type] = sample_list # add st type and its samples to a dictionary
        else:
            pass
    return st_dict

def create_sample_sheets(st_dict, samplesheet):
    """Create a samplesheet with the assemblies for each ST type."""
    for key, sample_list in st_dict.items():
        with open(key +"_samplesheet.csv", 'w') as new_samplesheet: # create a new sample sheet for each ST
            new_samplesheet.write('sample,seq_type,assembly\n') #write the header
            for sample in sample_list: # for each sample that is part of the ST
                with open(samplesheet, 'r') as f: # read the orginal griphin samplesheet
                    for line in f:
                        if sample in line:
                            sample = line.split(',')[0]
                            assembly = line.split(',')[1].strip() + "/Assembly/" + sample + ".filtered.scaffolds.fa.gz" + "\n"
                            new_samplesheet.write(sample + "," + key + "," + assembly)

def main():
    args = parseArgs()
    st_dict = get_st_groups(args.griphin_report) # open the excel sheet get sample names organized by their ST types
    create_sample_sheets(st_dict, args.samplesheet) # go back to the samplesheet and keep only lines this matching samplenames

if __name__ == '__main__':
    main()
