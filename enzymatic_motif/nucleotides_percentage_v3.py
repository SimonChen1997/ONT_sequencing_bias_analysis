import csv
import re
import sys
import os
import itertools
import argparse

def main():
    parser = argparse.ArgumentParser(description='A Python tool for calculating '
                                                 'the percentage of bases in each position '
                                                 'The reads for using this method should be over 100bp')
    parser.add_argument('-p','--path', type=str, required=True,
                        help='a path where the fasta files are stored')
    parser.add_argument('-l', '--long', type=int, required=True,
                        help='a length of subreads')                     
    args = parser.parse_args()
    
    os.chdir(args.path)
    if args.path is not None and os.path.exists(args.path):
        filelist = getfile(args.path)
        for filename in filelist:
            base_percentages(filename, args.long)
    else:
        print("Error: Path argument is missing.")

def getfile(path):
    filelist = []
    for name in os.listdir(path): # extract all the fastq file name in the working directory
        if name.endswith('.fasta'):
            filelist.append(name)
    return filelist

def base_percentages(fasta_file,long):
    ## pre-create a base count dictionary
    base_counts = [{"A": 0, "C": 0, "G": 0, "T": 0, "N": 0} for position in range(long)]

    sequencing_kit = ""
    if "LSK" in fasta_file:
        sequencing_kit = "Ligation_kit"
    elif "RBK" in fasta_file:
        sequencing_kit = "Rapid_kit"

    # Read the FASTA file
    with open(fasta_file, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                # Process the previous sequence
                if sequence:
                    for i, base in enumerate(sequence[:long]):
                        base_counts[i][base] += 1
                    sequence = ""
            else:
                # Accumulate the sequence
                sequence += line.strip()

        # Process the last sequence
        for i, base in enumerate(sequence[:long]):
            base_counts[i][base] += 1

    # Calculate percentages for each position
    total_bases = [sum(base_counts[i].values()) for i in range(long)]
    base_percentages = [{base: "{:.2f}".format(count / total_bases[i] * 100, 4) for base, count in base_counts[i].items()} for i in range(long)]

    # Write base percentages to a CSV file
    csv_filename = fasta_file.replace(".fasta", "_base_percentages.csv")
    sample_id=fasta_file.replace(".fasta","")
    with open(csv_filename, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["sample_id"]+["sequencing_kit"] + ["bases"] + ["percentages"] + ["positions"])  # Column names
        for i in range(long):
            writer.writerow([sample_id] + [sequencing_kit] + ["A"] + [base_percentages[i]["A"]] + [i+1])  # A row
        for i in range(long):
            writer.writerow([sample_id] + [sequencing_kit] + ["T"] + [base_percentages[i]["T"]] + [i+1])  # T row
        for i in range(long):
            writer.writerow([sample_id] + [sequencing_kit] + ["C"] + [base_percentages[i]["C"]] + [i+1])  # C row
        for i in range(long):
            writer.writerow([sample_id] + [sequencing_kit] + ["G"] + [base_percentages[i]["G"]] + [i+1])  # G row
        for i in range(long):
            writer.writerow([sample_id] + [sequencing_kit] + ["N"] + [base_percentages[i]["N"]] + [i+1])  # N row

if __name__ == "__main__":
    main()
