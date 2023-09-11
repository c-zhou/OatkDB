#!/usr/bin/env python

import argparse
import gzip
import re
from collections import defaultdict
from io import StringIO
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

def find_anticodon(feature):
    anticodon = feature.qualifiers.get("anticodon", [""])[0]
    if anticodon and re.match(r'^[ACGU]{3}$', anticodon):
        return anticodon
    else:
        # try to find anticodon in /note
        note = feature.qualifiers.get("note", [""])[0]
        if note:
            match = re.match(r'.*anticodon:([ACGU]{3}).*', note)
            if match:
                # anticodon found
                return match.group(1)
    return ""

# A custom function to check if a character is valid for a gene name
# for most commonly used gene names
def is_valid_char(char):
    return char.isalnum() or char in ('_', '-', '.')

# Extract the start and end positions
def location_string(location):
    return f"{location.nofuzzy_start}-{location.nofuzzy_end}"

def extract_genes(record, genetic_code, include_tRNA, include_rRNA, check_codon, cds_translation):
    genes = []

    # update genetic code if translation table is recorded
    for feature in record.features:
        if feature.type == "CDS":
            transl_table = feature.qualifiers.get("transl_table", [""])[0]
            if transl_table:
                # use genetic code recorded in the genbank file
                genetic_code = int(transl_table)
                break
    codon_table = CodonTable.generic_by_id[genetic_code].forward_table

    for feature in record.features:
        if feature.type == "CDS" or (feature.type == "tRNA" and include_tRNA) or (feature.type == "rRNA" and include_rRNA):
            gene = feature.qualifiers.get("gene", [""])[0]
            try:
                sequence = str(feature.extract(record.seq))
            except Exception as e:
                sequence = ""

            if feature.type == "tRNA":
                if not re.match(r'^trn[a-z]{0,1}[ACDEFGHIKLMNPQRSTVWY]{1}-[ACGU]{3}$', gene):
                    # to rescue tRNA genes with informal names
                    if re.match(r'^trn[a-z]{0,1}[ACDEFGHIKLMNPQRSTVWY]{1}\([ACGU]{3}\)$', gene):
                        # change trn*(*) to trn*-*
                        gene = re.sub(r'\((.*?)\)', r'-\1', gene)
                    if re.match(r'^trn[a-z]{0,1}[ACDEFGHIKLMNPQRSTVWY]{1}_[ACGU]{3}$', gene):
                        # change trn*_* to trn*-*
                        gene = gene.replace("_", "-")
                    if re.match(r'^trn[a-z]{0,1}[ACDEFGHIKLMNPQRSTVWY]{1}$', gene):
                        # missing anticodon
                        anticodon = find_anticodon(feature)
                        if anticodon:
                            gene = gene + "-" + anticodon;

                match = re.match(r'^trn[a-z]{0,1}([ACDEFGHIKLMNPQRSTVWY]{1})-([ACGU]{3})$', gene)
                if not match:
                    # gene name cannot be formatted
                    gene = ""
                elif check_codon:
                    # check if in the codon table
                    aa = match.group(1)
                    anticodon = match.group(2)
                    if codon_table.get(str(Seq(anticodon).reverse_complement())) != aa:
                        # anticodon and AA do not match
                        gene = ""

            if feature.type == "rRNA":
                # fine to change names to lowercase for rRNA genes
                gene = gene.lower()
                if gene.endswith("s") or gene.endswith("l"):
                    # change rrn[4.5|5|16|23][s|S] and rrn[4.5|5|16|23][l|L] to rrn[4.5|5|16|23]
                    # for small subunit and large subunit
                    gene = gene[:-1]
                match = re.match(r'^(\d+(\.\d+)?)rrn$', gene)
                if match:
                    # change [4.5|5|16|23]rrn to rrn[4.5|5|16|23]
                    gene = "rrn" + match.group(1)
                if re.match(r'^\d+(\.\d+)?$', gene):
                    # change [4.5|5|16|23] to rrn[4.5|5|16|23]
                    gene = "rrn" + gene
                if not re.match(r'^rrn\d+(\.\d+)?$', gene):
                    # gene name cannot be formatted
                    gene = ""

            if feature.type == "CDS" and cds_translation:
                translation = feature.qualifiers.get("translation", [""])[0]
                
                # If translation is not available, calculate it using the specified genetic code
                if not translation:
                    seq_obj = Seq(sequence)
                    try:
                        translation = str(seq_obj.translate(table=genetic_code, cds=True))
                    except CodonTable.TranslationError:
                        translation = ""

                sequence = translation
            
            if gene and all(is_valid_char(char) for char in gene) and sequence:
                genes.append([feature.type, gene, record.id, location_string(feature.location), sequence])
    
    return genes

def main(input_genbank_file, output_file, summary_file, genetic_code, no_rRNA, no_tRNA, check_codon, translate):
    gene_list = []
    gene_counts = defaultdict(int)
    '''
    for record in SeqIO.parse(input_genbank_file, "genbank"):
        genes = extract_genes(record, genetic_code, not no_tRNA, not no_rRNA, check_codon, translate)
        gene_list.extend(genes)
    '''
    # this is to replace the naive SeqIO parse to deal with GenBank file parsing errors
    with open(input_genbank_file, "r") as genbank_file:
        gb_file_buff = []
        for line in genbank_file:
            line = line.rstrip()
            if line:
                gb_file_buff.append(line)
                if line == "//":
                    # Genbank record end reached
                    # pack record lines into a StringIO object for SeqIO parsing
                    string_io = StringIO("\n".join(gb_file_buff))
                    try:
                        record = SeqIO.read(string_io, "genbank")
                        genes = extract_genes(record, genetic_code, not no_tRNA, not no_rRNA, check_codon, translate)
                        gene_list.extend(genes)
                    except Exception as e:
                        print(f"Error parsing record: {gb_file_buff[0]}")
                    string_io.close()
                    gb_file_buff.clear()

    # find the most common name for each case-insensitive gene name
    name_counts = {}
    for _, gene_name, _, _, _ in gene_list:
        name_lower = gene_name.lower()
        if name_lower not in name_counts:
            name_counts[name_lower] = defaultdict(int)
        name_counts[name_lower][gene_name] += 1

    name_dict = defaultdict(str)
    for name_lower, name_count in name_counts.items():
        name_dict[name_lower] = max(name_count, key=name_count.get)
    
    # change gene name to the most commonly used one
    for gene in gene_list:
        gene[1] = name_dict[gene[1].lower()]

    # Sort the genes by gene type and then by gene name
    sorted_genes = sorted(gene_list, key=lambda x: (x[0], x[1], x[2]))

    # Check if the output file name ends with ".gz" and write a gzipped file if it does
    if output_file.endswith(".gz"):
        with gzip.open(output_file, "wt") as output:
            output.write("#Type\tGene\tSequence\tLocation\tSequence\n")
            for gene_type, gene_name, sequence_id, location, sequence in sorted_genes:
                output.write(f"{gene_type}\t{gene_name}\t{sequence_id}\t{location}\t{sequence}\n")
    else:
        with open(output_file, "w") as output:
            output.write("#Type\tGene\tSequence\tLocation\tSequence\tTranslation\n")
            for gene_type, gene_name, sequence_id, location, sequence in sorted_genes:
                output.write(f"{gene_type}\t{gene_name}\t{sequence_id}\t{location}\t{sequence}\n")

    # Count genes by type and name for the summary file
    for gene_type, gene_name, _, _, _ in sorted_genes:
        gene_counts[f"{gene_type}_{gene_name}"] += 1

    # Sort the genes by count for the summary file
    sorted_gene_counts = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)

    # Write the summary file, sorted by counts
    with open(summary_file, "w") as summary_output:
        summary_output.write("#Type\tGene\tCount\n")
        for gene, count in sorted_gene_counts:
            tgene=gene.replace('_', '\t', 1)
            summary_output.write(f"{tgene}\t{count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract coding genes (CDS), tRNA, and rRNA genes from a GenBank file with multiple sequences, store them in a list, sort by gene types and names, and output in tabular format with a summary file sorted by counts. Optionally, specify a genetic code table and decide whether to include tRNA and rRNA, and to write translation sequences for coding genes.")
    parser.add_argument("input_genbank_file", help="Path to the input GenBank file")
    parser.add_argument("output_file", help="Path to the output file in tabular format (can end with '.gz' for gzipped output)")
    parser.add_argument("summary_file", help="Path to the summary file for gene counts (sorted by counts)")
    parser.add_argument("--genetic_code", type=int, default=1, help="Genetic code table number for translation (default: 1). For details on available code tables, visit https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
    parser.add_argument("--no_rRNA", action="store_true", help="Exclude rRNA genes")
    parser.add_argument("--no_tRNA", action="store_true", help="Exclude tRNA genes")
    parser.add_argument("--check_codon", action="store_true", help="Check codon consistency for tRNA genes")
    parser.add_argument("--translate", action="store_true", help="Write translation sequences for CDS")

    args = parser.parse_args()
    main(args.input_genbank_file, args.output_file, args.summary_file, args.genetic_code, args.no_rRNA, args.no_tRNA, args.check_codon, args.translate)

