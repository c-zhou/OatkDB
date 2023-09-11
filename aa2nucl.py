#!/usr/bin/env python

import argparse
import math
import sys
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio import SeqIO

def calculate_lcm_of_list(numbers):
    # Initialize the LCM to the first number in the list
    lcm = numbers[0]
    # Iterate through the list and calculate the LCM
    for num in numbers[1:]:
        # Calculate the LCM using the formula: LCM(a, b) = (a * b) / GCD(a, b)
        lcm = (lcm * num) // math.gcd(lcm, num)
    return lcm

class Transl:
    # amino acid to code table used for conversion
    def __init__(self, codon_table_id=1):
        self.codon_table_id = codon_table_id
        # LCM or sequence number
        self.seq_number = 0
        # Initialize an empty dictionary to store the mapping
        self.amino_acid_codon_dict = {'*' : ['---']}
        self.initialize()

    def initialize(self):
        # Get the specified codon table by ID
        codon_table = CodonTable.unambiguous_dna_by_id[self.codon_table_id]
        # Iterate through the codon table's forward table (codon, amino_acid)
        for codon, amino_acid in codon_table.forward_table.items():
            if amino_acid in self.amino_acid_codon_dict:
                self.amino_acid_codon_dict[amino_acid].append(codon)
            else:
                self.amino_acid_codon_dict[amino_acid] = [codon]
        # Collect dictionary size for each amino acid
        # and calculate the LCM
        lengths_list = [len(value) for value in self.amino_acid_codon_dict.values()]
        self.seq_number = calculate_lcm_of_list(lengths_list)
        # Expand the codon list for each amino acid
        for amino_acid in self.amino_acid_codon_dict:
            n = int(self.seq_number / len(self.amino_acid_codon_dict[amino_acid]))
            self.amino_acid_codon_dict[amino_acid] *= n

    def protein_to_nucleotide(self, protein_sequence):
        try:
            # Create an empty nucleotide sequence
            seq_length = len(protein_sequence)
            codon_mat = [[None] * seq_length for _ in range(self.seq_number)]
            # Iterate through the amino acids in the protein sequence
            seq_pos = 0
            for amino_acid in protein_sequence:
                if amino_acid not in self.amino_acid_codon_dict:
                    amino_acid = '*'
                codons = self.amino_acid_codon_dict[amino_acid]
                for s in range(self.seq_number):
                    codon_mat[s][seq_pos] = codons[s]
                seq_pos += 1
            nucleotide_sequences = [''.join(codon_row) for codon_row in codon_mat]
            return nucleotide_sequences
        except Exception as e:
            print(f"Error: {str(e)}")
            return None
    
    def write_protein_to_nucleotide(self, protein_sequence, output_handle, seq_id):
        try:
            # Create an empty nucleotide sequence
            seq_length = len(protein_sequence)
            # Iterate through the amino acids in the protein sequence
            for i in range(self.seq_number):
                output_handle.write(f">{seq_id}-{i+1}\n")
                for j in range(seq_length):
                    amino_acid = protein_sequence[j]
                    if amino_acid not in self.amino_acid_codon_dict:
                        amino_acid = '*'
                    output_handle.write(self.amino_acid_codon_dict[amino_acid][i])
                output_handle.write("\n")
        except Exception as e:
            print(f"Error: {str(e)}")

def convert_protein_sequences(input_file, output_file, codon_table_id=1):
    try:
        # Make translate table
        trans = Transl(codon_table_id)
        # Open the input and output files
        if input_file and input_file != "-":
            input_handle = open(input_file, "r")
        else:
            input_handle = sys.stdin
        if output_file and output_file != "-":
            output_handle = open(output_file, 'w')
        else:
            output_handle = sys.stdout
        
        for record in SeqIO.parse(input_handle, "fasta"):
            protein_sequence = record.seq.upper()
            nucleotide_sequences = trans.protein_to_nucleotide(protein_sequence)
            if nucleotide_sequences:
                for s in range(trans.seq_number):
                    # Write the nucleotide sequence to the output file
                    output_handle.write(f">{record.id}-{s+1}\n")
                    output_handle.write(f"{nucleotide_sequences[s]}\n")
            #trans.write_protein_to_nucleotide(protein_sequence, output_handle, record.id)

        # Close the input and output file handles
        if input_file:
            input_handle.close()
        if output_file:
            output_handle.close()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Convert protein sequence alignments to nucleotide sequence alignments for HMM profile construciton. See Fischer, Carlos N., et al 2015.")
    parser.add_argument("input_file", nargs="?", default=None, help="Input FASTA file with protein sequences (default: stdin)")
    parser.add_argument("output_file", nargs="?", default=None, help="Output FASTA file to write nucleotide sequences (default: stdout)")
    parser.add_argument("--codon_table_id", type=int, default=1, help="Codon table ID (default: 1)")
    
    args = parser.parse_args()
    
    convert_protein_sequences(args.input_file, args.output_file, args.codon_table_id)

