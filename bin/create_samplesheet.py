#!/usr/bin/env python3

import sys
import glob
import os
import argparse
import csv

##Given a samplesheet from GRiPHin's wf a new samplesheet
##Usage: >python create_samplesheet.py -s ./samplesheet.csv -a ../PHX/phoenix/assets/databases/ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a samplesheet with sample,directory columns.')
    parser.add_argument('-d', '--directory', default=None, required=False, dest='directory', help='Will create a samplesheet for all samples in the directory (expects PHoeNIx style directory structure).')
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def create_samplesheet(directory):
    """Function will create a samplesheet from samples in a directory if -d argument passed."""
    directory = os.path.abspath(directory) # make sure we have an absolute path to start with
    with open("Directory_samplesheet.csv", "w") as samplesheet:
        samplesheet.write('sample,directory\n')
        dirs = os.listdir(directory)
        # If there are any new files added to the top directory they will need to be added here or you will get an error
        skip_list_a = glob.glob(directory + "/*_GRiPHin_Summary.*") # for if griphin is run on a folder that already has a report in it
        skip_list_a = [ gene.split('/')[-1] for gene in skip_list_a ]  # just get the excel name not the full path
        skip_list_b = ["BiosampleAttributes_Microbe.1.0.xlsx", "Sra_Microbe.1.0.xlsx", "Phoenix_Summary.tsv", "pipeline_info", "GRiPHin_Summary.xlsx", "GRiPHin_Summary.tsv", "multiqc", "samplesheet_converted.csv", "Directory_samplesheet.csv", "sra_samplesheet.csv"]
        skip_list = skip_list_a + skip_list_b
        for sample in dirs:
            if sample not in skip_list:
                #with open("Directory_samplesheet.csv", "a") as samplesheet:
                    if directory[-1] != "/": # if directory doesn't have trailing / add one
                        directory = directory + "/"
                    #print(sample + "," + directory + sample + '\n')
                    samplesheet.write(sample + "," + directory + sample + '\n')
                    #print(directory)
    samplesheet = "Directory_samplesheet.csv"
    return samplesheet

def main():
    args = parseArgs()
    # If a directory is given then create a samplesheet from it if not use the samplesheet passed
    if args.directory !=None:
        samplesheet = create_samplesheet(args.directory)
    else:
        sys.exit(CRED + "You MUST pass a stop directory of PHoeNIx output to create a samlesheet.\n" + CEND)

if __name__ == '__main__':
    main()