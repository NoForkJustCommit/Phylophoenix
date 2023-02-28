#!/usr/bin/env python3

import pandas as pd
import argparse
import csv
import glob

##Makes a new samplesheet for assemblies based on their st type. Expects that assemblies are in an output PHoeNIx folder
##Usage: >python get_st_types.py -s GRiPHin_samplesheet.csv -g GRiPHin_Report.xlsx
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--samplesheet', default=None, required=False, dest='samplesheet', help='GRiPHin samplesheet of sample,directory in csv format.')
    return parser.parse_args()

def create_sample_sheets(samplesheet):
    """Create a samplesheet with the assemblies for each Seq Type. Also, creates samplesheet to run SNVPhyl for each Seq Type."""
    complete_list = []
    seq_type = "All_STs"
    with open("SNVPhyl_" + seq_type +"_samplesheet.csv", "a") as st_snv_samplesheet: # create a new sample sheet for each ST that can be used by snvphyl
        st_snv_samplesheet.write('id,directory') #write the header
    df = pd.read_csv(samplesheet, sep=',', header=0)
    sample_list = df["sample"].tolist()
    for sample in sample_list: # for each sample that is part of the ST
        with open(samplesheet, 'r') as f: # read the orginal griphin samplesheet
            for line in f:
                if str(sample) in line:
                    with open("SNVPhyl_" + seq_type +"_samplesheet.csv", "a") as st_snv_samplesheet: # this create a file with headers
                        st_snv_samplesheet.write("\n" + line.strip('\n'))
                    assembly = line.split(',')[1].strip() + "/Assembly/" + str(sample) + ".filtered.scaffolds.fa.gz"
                    complete_list.append(assembly)
    with open(seq_type +"_samplesheet.csv", 'w') as new_samplesheet: # create a new sample sheet for each ST
        new_samplesheet.write("id,seq_type,assembly_1,assembly_2") #write the header
        combinations = [(a, b) for idx, a in enumerate(complete_list) for b in complete_list[idx + 1:]]
        for combo in combinations:
            sample1 = combo[0].split('/')[-1].replace(".filtered.scaffolds.fa.gz","")
            sample2 = combo[1].split('/')[-1].replace(".filtered.scaffolds.fa.gz","")
            new_samplesheet.write("\n" + sample1 + "_" + sample2 + "," + seq_type + "," + str(combo[0]) + "," + str(combo[1]))

def main():
    args = parseArgs()
    create_sample_sheets(args.samplesheet) # go back to the samplesheet and keep only lines this matching samplenames

if __name__ == '__main__':
    main()
