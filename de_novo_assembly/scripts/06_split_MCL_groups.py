import os
import sys
from Bio import SeqIO

# Usage: python script.py input_ids.txt input_fasta.fa
input_ids_file = sys.argv[1]
input_fasta_file = sys.argv[2]

output_directory = "MCL_sequence_groups"
os.makedirs(output_directory, exist_ok=True)

# Read fasta file and create a dictionary with sequence IDs as keys
records_dict = SeqIO.to_dict(SeqIO.parse(input_fasta_file, "fasta"))

# Process the group file
with open(input_ids_file, 'r') as group_file:
    for index, line in enumerate(group_file):
        group_ids = line.strip().split('\t')
        group_sequences = []

        # Retrieve sequences by their IDs
        for seq_id in group_ids:
            if seq_id in records_dict:
                group_sequences.append(records_dict[seq_id])

        # Write sequences to a new fasta file
        output_file = os.path.join(output_directory, f"Group_{index+1}.fasta")
        with open(output_file, 'w') as outfile:
            SeqIO.write(group_sequences, outfile, "fasta")
