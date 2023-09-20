#!/usr/bin/env python

import argparse
import gzip
import sys
import re
from collections import defaultdict
from io import StringIO
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

class Alias:
    # alias gene name conversion

    def __init__(self, target_string, regex_pattern):
        self.target_string = target_string
        self.regex_pattern = re.compile(regex_pattern)
        self.check_compatibility()

    def check_compatibility(self):
        max_group_index = max((int(match) for match in re.findall(r'\((\d+)\)', self.target_string)), default = 0);
        if max_group_index > self.regex_pattern.groups:
            print(f"Error: undefined caputre group in target string: '{self.target_string}'")
            print(f"       regex pattern: '{self.regex_pattern.pattern}'")
            sys.exit(1)

    def alias_string_format(self, input_string):
        # Use re.search() to find the pattern within the target_string
        '''
        Examples:
        input_string    target_string   regex_pattern (pre-compiled)        output_string
        COI             COX1            ^[Cc][Oo][Xx]?(1|I)$                COX1
        COX2            COX2            ^[Cc][Oo][Xx]?(2|II)$               COX2
        cox3            COX3            ^[Cc][Oo][Xx]?(3|III)$              COX3
        Cyt-b           cytb            ^[Cc][Yy][Tt]?o?-?[B|b]$            cytb
        ndh4L           ND4L            ^[Nn][Aa]?[Dd][Hh]?4[Ll]$           ND4L
        ND5             ND(1)           ^[Nn][Aa]?[Dd][Hh]?(\d{1})$         ND5
        ndhB            ndh(1)          ^[Nn][Aa]?[Dd][Hh]?([a-zA-Z]{1})$   ndhB

        Example usage:
        alias_pattern = Alias("cytb", "^[Cc][Yy][Tt]?o?-?[B|b]$")
        result = alias_pattern.alias_string_format("Cytb")
        '''
        # Use re.search() to find the pattern within the input_string
        match = self.regex_pattern.match(input_string)

        if match:
            # Extract all matched groups
            matched_groups = match.groups()
            target_string = str(self.target_string)
            # Replace "(1)", "(2)", etc., in the target_string with the matched groups
            for i, group in enumerate(matched_groups, 1):
                target_string = target_string.replace(f'({i})', group)

            return True, target_string
        else:
            # Return the original input_string if no match is found
            return False, input_string

def find_anticodon(feature):
    anticodon = feature.qualifiers.get("anticodon", [""])[0]
    if anticodon and re.match(r'^[ACGUacgu]{3}$', anticodon):
        return anticodon.upper()
    else:
        # try to find anticodon in /note
        note = feature.qualifiers.get("note", [""])[0]
        if note:
            match = re.match(r'.*anticodon:([ACGUacgu]{3}).*', note)
            if match:
                # anticodon found
                return match.group(1).upper()
    return ""

# A custom function to check if a character is valid for a gene name
# for most commonly used gene names
def is_valid_char(char):
    return char.isalnum() or char in ('_', '-', '.')

# Extract the start and end positions
def location_string(location):
    return f"{location.nofuzzy_start}-{location.nofuzzy_end}"

def extract_genes(record, genetic_code, include_tRNA, include_rRNA, include_ORF, check_codon, cds_translation):
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
                matched = False
                aminoacid = None
                anticodon = None

                # match trnM-GGG trnM_GGG trnM(GGG) TrNM[GGG] trnM trnfM etc
                match =  re.match(r'^([t|T][r|R][n|N])([a-z]{0,1})([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]{1})[-|_]?[\(|\]]?([ACGUacgu]{3})?[\)|\]]?$', gene)
                if match:
                    aminoacid = match.group(3).upper()
                    anticodon = match.group(4)
                    # to rescue tRNA genes without anticodon
                    if not anticodon:
                        anticodon = find_anticodon(feature)
                    else:
                        anticodon = anticodon.upper()
                    if anticodon:
                        matched = True
                        gene = f"{match.group(1).lower()}{match.group(2)}{aminoacid}-{anticodon}"

                if not matched:
                    # gene name cannot be formatted
                    gene = ""
                elif check_codon:
                    # check if in the codon table
                    if codon_table.get(str(Seq(anticodon).reverse_complement())) != aminoacid:
                        # anticodon and AA do not match
                        gene = ""

            if feature.type == "rRNA":
                # fine to change names to lowercase for rRNA genes
                gene = gene.lower()
                # match [4.5|5|16|23][s?]rrn
                match = re.match(r'^(\d+(\.\d+)?)s?rrn$', gene)
                if not match:
                    # match [4.5|5|16|23][s?]
                    match = re.match(r'^(\d+(\.\d+)?)s?$', gene)
                if not match:
                    # match rrn[4.5|5|16|23][s?]
                    match = re.match(r'^rrn(\d+(\.\d+)?)s?$', gene)
                if not match:
                    # gene name cannot be formatted
                    # to rescue rRNA genes with product information
                    product = feature.qualifiers.get("product", [""])[0]
                    if product:
                        product = product.lower()
                        match = re.search(r'(\d+(\.\d+)?)s?[\s\t]*ribosom(al|e)[\s\t]*rna', product)
                    # to rescue rRNA genes with note information
                    if not match:
                        note = feature.qualifiers.get("note", [""])[0]
                        if note:
                            note = note.lower()
                        match = re.search(r'(\d+(\.\d+)?)s?[\s\t]*ribosom(al|e)[\s\t]*rna', note)

                if match:
                    gene = f"rrn{match.group(1)}"
                else:
                    gene = ""

            if feature.type == "CDS":
                # TODO is this safe for ORF checking?
                if not include_ORF and re.match(r'.*[O|o][R|r][F|f].*', gene):
                    gene = ""

                if gene and cds_translation:
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

def parse_alias_file(input_alias_file):
    aliases = []
    with open(input_alias_file, "r") as alias_file:
        for line in alias_file:
            line = line.strip()

            # Skip empty lines and lines that start with "#"
            if not line or line.startswith("#"):
                continue

            columns = re.split(r'\s+', line.strip())
            # Ensure there are at least two columns before extracting
            if len(columns) < 2:
                print(f"Error: not a valild alias rule: '{self.target_string}'")
                print(f"       at lease two columns required")
                sys.exit(1)
            
            # Append the new Alias to the data list
            aliases.append(Alias(columns[0], columns[1]))
    
    return aliases

def parse_genbank_file(input_genbank_file, output_file, summary_file, genetic_code, gene_alias, no_rRNA, no_tRNA, include_ORF, check_codon, translate):
    gene_list = []
    '''
    for record in SeqIO.parse(input_genbank_file, "genbank"):
        genes = extract_genes(record, genetic_code, not no_tRNA, not no_rRNA, include_ORF, check_codon, translate)
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
                        genes = extract_genes(record, genetic_code, not no_tRNA, not no_rRNA, include_ORF, check_codon, translate)
                        gene_list.extend(genes)
                    except Exception as e:
                        print(f"Error parsing record: {gb_file_buff[0]}")
                    string_io.close()
                    gb_file_buff.clear()
    
    # read gene name aliases
    aliases = parse_alias_file(gene_alias)
    
    # convert gene name to aliases
    name_aliases = defaultdict(str)
    gene_set = set(map(lambda item: item[1], gene_list))
    for gene in gene_set:
        for alias in aliases:
            matched, name_alias = alias.alias_string_format(gene)
            if matched:
                break
        name_aliases[gene] = name_alias

    # find the most common name for each case-insensitive gene alias name
    name_counts = {}
    for _, gene_name, _, _, _ in gene_list:
        name_lower = name_aliases[gene_name].lower()
        if name_lower not in name_counts:
            name_counts[name_lower] = defaultdict(int)
        name_counts[name_lower][gene_name] += 1

    name_dict = defaultdict(str)
    for name_lower, name_count in name_counts.items():
        name_dict[name_lower] = max(name_count, key=name_count.get)
    
    # change gene name to the most commonly used one
    for gene in gene_list:
        gene[1] = name_dict[name_aliases[gene[1]].lower()]

    # Sort the genes by gene type and then by gene name
    sorted_genes = sorted(gene_list, key=lambda x: (x[0], x[1], x[2]))
    
    # Update sequence gene counts
    seq_gene_counts = defaultdict(int)
    # Count genes by sequence
    for _, _, sequence_id, _, _ in sorted_genes:
        seq_gene_counts[sequence_id] += 1

    # Check if the output file name ends with ".gz" and write a gzipped file if it does
    if output_file.endswith(".gz"):
        with gzip.open(output_file, "wt") as output:
            output.write("#Type\tGene\tAccession\tCount\tLocation\tSequence\n")
            for gene_type, gene_name, sequence_id, location, sequence in sorted_genes:
                output.write(f"{gene_type}\t{gene_name}\t{sequence_id}\t{seq_gene_counts[sequence_id]}\t{location}\t{sequence}\n")
    else:
        with open(output_file, "w") as output:
            output.write("#Type\tGene\tAccession\tCount\tLocation\tSequence\n")
            for gene_type, gene_name, sequence_id, location, sequence in sorted_genes:
                output.write(f"{gene_type}\t{gene_name}\t{sequence_id}\t{seq_gene_counts[sequence_id]}\t{location}\t{sequence}\n")

    gene_counts = defaultdict(int)
    # Count genes by type and name for the summary file
    for gene_type, gene_name, _, _, _ in sorted_genes:
        gene_counts[f"{gene_type}_{gene_name}"] += 1

    # Sort the genes by count for the summary file
    sorted_gene_counts = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)

    # Write the summary file, sorted by counts
    with open(summary_file, "w") as summary_output:
        summary_output.write("#Type\tGene\tName\tCount\n")
        for gene, count in sorted_gene_counts:
            tgene=gene.split('_', 1)
            summary_output.write(f"{tgene[0]}\t{tgene[1]}\t{tgene[1]}\t{count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract coding genes (CDS), tRNA, and rRNA genes from a GenBank file with multiple sequences, store them in a list, sort by gene types and names, and output in tabular format with a summary file sorted by counts. Optionally, specify a genetic code table and decide whether to include tRNA and rRNA, and to write translation sequences for coding genes.")
    parser.add_argument("input_genbank_file", help="Path to the input GenBank file")
    parser.add_argument("output_file", help="Path to the output file in tabular format (can end with '.gz' for gzipped output)")
    parser.add_argument("summary_file", help="Path to the summary file for gene counts (sorted by counts)")
    parser.add_argument("--genetic_code", type=int, default=1, help="Genetic code table number for translation (default: 1). For details on available code tables, visit https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
    parser.add_argument("--gene_alias", type=str, default="", help="Gene name aliases file")
    parser.add_argument("--no_rRNA", action="store_true", help="Exclude rRNA genes")
    parser.add_argument("--no_tRNA", action="store_true", help="Exclude tRNA genes")
    parser.add_argument("--include_ORF", action="store_true", help="Include open reading frames")
    parser.add_argument("--check_codon", action="store_true", help="Check codon consistency for tRNA genes")
    parser.add_argument("--translate", action="store_true", help="Write translation sequences for CDS")
    
    args = parser.parse_args()
    parse_genbank_file(args.input_genbank_file, args.output_file, args.summary_file, args.genetic_code, args.gene_alias, args.no_rRNA, args.no_tRNA, args.include_ORF, args.check_codon, args.translate)

