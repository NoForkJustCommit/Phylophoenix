#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import csv
import glob
import re

##Makes a new samplesheet for assemblies based on their st type. Expects that assemblies are in an output PHoeNIx folder
##Usage: >python get_reference_seq.py -s GRiPHin_samplesheet.csv -i sample1_sample2.tsv -t seq_type
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to return the full file path for the centroid assembly.')
    parser.add_argument('-i', '--input', required=True, dest='input', help='txt file that has mash distances on each line.')
    parser.add_argument('-s', '--samplesheet', required=True, dest='samplesheet', help='GRiPHin_samplesheet.csv.')
    parser.add_argument('-t', '--seq_type', required=True, dest='seq_type', help='sequence type.')
    return parser.parse_args()

def average_mash(input_mash_list):
    """Takes in a cat list of mash distances and then makes a list with the average value for each entry"""
    #Read in tsv file as a pandas datafame
    mash_df = pd.read_csv(input_mash_list, sep='\t', names=["Reference_ID", "Query_ID", "Mash_distance", "P_value", "Matching_hashes"], index_col=False, dtype={'Mash_distance': np.float64, 'P_value': np.float64})
    #clean up sample names
    mash_df['Reference_ID'] = mash_df['Reference_ID'].str.replace(".filtered.scaffolds.fa.gz", "", regex=False)
    mash_df['Query_ID'] = mash_df['Query_ID'].str.replace(".filtered.scaffolds.fa.gz", "", regex=False)
    #Get unique sample names
    samples = mash_df['Query_ID'].unique()
    #filter dataframe to only have one sample
    count = 0
    for sample in samples:
        filtered_df = mash_df.loc[(mash_df['Reference_ID'] == sample) | (mash_df['Query_ID'] == sample)]
        ave_dis = filtered_df['Mash_distance'].mean() #get average of distances
        if count == 0: # initial setting of low_mean and sample
            low_mean = ave_dis
            low_mean_sample = sample
            count = count + 1 # turn this off
        elif ave_dis < low_mean: # set new low mean if the ave_dis for this sample is lower
            low_mean = ave_dis
            low_mean_sample = sample
        elif ave_dis == low_mean:
            print("Two samples have the same mean distance.")
    return low_mean_sample

def get_assembly_path(seq_type,samplesheet,low_mean_sample):
    with open(samplesheet, "r") as sheet:
        header = next(sheet) #skip header
        for line in sheet:
            match = re.search(low_mean_sample, line)
            if match is not None:
                print(low_mean_sample)
                print(line.split(",")[1])
                path_to_centroid_dir = line.split(",")[1].strip("\n")
                if path_to_centroid_dir.endswith("/"):
                    path_to_centroid_dir.strip("/")
                path_to_centroid = path_to_centroid_dir + "/assembly/" + match[0] + ".filtered.scaffolds.fa.gz"
    with open("path_to_" + seq_type + "_centroid.csv", "w") as output:
        output.write(path_to_centroid)

def main():
    args = parseArgs()
    low_mean_sample = average_mash(args.input) # open the excel sheet get sample names organized by their ST types
    get_assembly_path(args.seq_type,args.samplesheet,low_mean_sample)

if __name__ == '__main__':
    main()