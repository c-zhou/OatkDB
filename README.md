# Oatk HMM profile database

## Overview
This repository hosts all versions of the HMM profile database used by the de novo organelle genome assembly tool [Oatk](https://github.com/c-zhou/oatk) and scripts for building these databases. The database files are inputs of [HMMER](http://hmmer.org/) for organelle sequence annotation. The initial version database ([v20230210](https://github.com/c-zhou/OatkDB/tree/main/v20230210)) was constructed using the pipeline described in [fppa](https://github.com/tolkit/fppa) and [fpma](https://github.com/tolkit/fpma). The scripts hosted in this repo aim to provide a more automated and user-friendly toolset for database construction.

## Installation
Download the source code from this repo or with `git clone https://github.com/c-zhou/OatkDB.git`

The tools are Shell scripts with some auxiliary Python scripts. You need a system supporting Shell/Bash script to run the tool. You also need [Python](https://www.python.org/) installed. The major dependency of the Python scripts is [Biopython](https://biopython.org/).

Other dependencies include:

* [seqtk](https://github.com/lh3/seqtk)
* [mafft](https://mafft.cbrc.jp/alignment/software/)
* [hmmer](http://hmmer.org/)
* [parallel](https://www.gnu.org/software/parallel/)
* [Entrez Direct: E-utilities on the Unix Command Line](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

These tools must be available in the system path, i.e., add them to the `$PATH` environmental variable. The program will check the availabilities of the dependencies before doing any actual work.

## Oatkdb
This is the major tool for building a database. It takes two positional parameters: the NCBI taxonomy ID of the species and the organelle type (mitochondria, chloroplast, or plastid). Here is how it works: (1) download all sequences in GenBank format from the NCBI database given the user-specified NCBI taxonomy and the organelle type; (2) parse the GenBank files of the source reference sequences to find core genes including protein-coding, rRNA and tRNA (mainly by the frequency seen in the source sequences) and extract corresponding sequences; (3) sequence clean to remove outlier/error-prone sequences considering sequence length, content, and divergence; (4) sequence multiple alignment, and finally (5) profile HMM construction.

It is worth noting that the whole pipeline is quite computationally intensive, especially the multiple alignment step, so it is critical to give as many CPUs as possible. This is controlled by the `-j` and `-t` options. Also, the tool can resume a failed run (`--resume` option), for which the intermediate/temporary files are needed, so it is highly recommended to keep them during the run (i.e., run the tool WITHOUT the `--clean` option). The intermediate files also contain many useful information. The temporary folder is also a mandatory input for `mancdb` (see below).

Here is the full list of the options, which can be shown with the command `oatkdb -h`.

```
  Program: oatkdb (build gene HMM profile database for an NCBI taxonomy)
  Version: 1.0

  Usage: oatkdb [options] ${taxid} [mitochondrion|chloroplast|plastid]
  Optional:
      -j|--jobs     INT    Number of parallel jobs to run (default 4).
      -t|--threads  INT    Number of threads to use for each parallel job (default 2).
      -p|--protein         Use protein instead of nucleotide sequence alignment for coding genes.
      -c|--codon    INT    NCBI genetic code table for translation (default 1).
      -T|--tmpdir   STR    Temp file directory (default auto).
      -f|--force           Overwrite existing files.
      -a|--alias    STR    Gene name alias file.
      -o|--output   STR    Output HMM profile database file name prefix (default auto).
         --max-ref  INT    Maximum number of NCBI reference sequences to download (default INF).
         --min-seq  INT    Minimum number of sequences to keep a gene (default 5).
         --max-seq  INT    Maximum number of sequences to use for HMM (default 10000).
         --no-trna         Do not build HMM profiles for tRNA genes.
         --no-rrna         Do not build HMM profiles for rRNA genes.
         --incl-orf        Include open reading frames.
         --resume          Resume a previously unfinished run.
         --clean           Clean temporary files when finished.
         --test            Run program in test mode.
         --log      STR    Log file (default stdout).
      -h|--help            Print this help message.
      -v|--version         Print version number.

  1. Refer to NCBI Taxonomy https://www.ncbi.nlm.nih.gov/taxonomy for taxid
  2. Refer to webpage https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for code table selection
     Commonly used genetic code tables include:
       [1]  The Standard Code
       [2]  The Vertebrate Mitochondrial Code
       [5]  The Invertebrate Mitochondrial Code
       [11] The Bacterial, Archaeal and Plant Plastid Code

  Example: oatkdb -j 4 -t 8 -c 11 -o angiosperms_pltd_v20230911 3398 chloroplast
```

## Seqdb
This tool is somewhat a subroutine of the `oatkdb` tool. It does part of the job that `oatkdb` does. It takes one positional parameter, i.e., the sequence file (in FASTA format) for a gene, and builds a profile HMM from the sequences. The work logic is similar to `oatkdb`: clean sequences, do multiple sequence alignment and build an HMM profile.

Here is the full list of the options, which can be shown with the command `seqdb -h`.

```
  Program: seqdb (build profile HMM for a gene from sequence file)
  Version: 1.0

  Usage: seqdb [options] sequence
  Optional:
      -n|--name     STR    Gene name to write to the HMM file (default auto).
      -t|--threads  INT    Number of threads to use for each parallel job (default 8).
      -p|--protein         The input is protein sequences.
      -c|--codon    INT    NCBI genetic code table for translation used for '-p' (default 1).
      -T|--tmpdir   STR    Temp file directory (default auto).
      -f|--force           Overwrite existing files.
      -o|--output   STR    Output HMM profile database file (default auto).
         --min-seq  INT    Minimum number of sequences to keep a gene (default 5).
         --max-seq  INT    Maximum number of sequences to use for HMM (default 10000).
         --resume          Resume a previously unfinished run.
         --clean           Clean temporary files when finished.
         --log      STR    Log file (default stdout).
      -h|--help            Print this help message.
      -v|--version         Print version number.

  1. Refer to NCBI Taxonomy https://www.ncbi.nlm.nih.gov/taxonomy for taxid
  2. Refer to webpage https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for code table selection
     Commonly used genetic code tables include:
       [1]  The Standard Code
       [2]  The Vertebrate Mitochondrial Code
       [5]  The Invertebrate Mitochondrial Code
       [11] The Bacterial, Archaeal and Plant Plastid Code

  Example: seqdb -t 12 -c 11 -n rpoc2 -o ./rpoc2.hmm rpo2.fna
```

## Mancdb
This tool is used for the final manual curation of the database. Even though `oatkdb` endeavours to build a clean core gene set, it is sometimes not perfect. `mancdb` provides an interface to manually curate (such as add and delete genes, and change gene names) the core gene set to be included in the database and rebuild the database. It takes three positional parameters: a gene list file where each row represents a gene to be included in the database, the temporary folder (used to run `oatkdb`) and the output database file name. The gene list file follows the same format as the intermediate file `DB.summary` which can be found in the `oatkdb` temporary folder. For each row, the first column (`Type`) and the second column (`Gene`) are used to located the corresponding sequence files for the gene and should be directly copied from `DB.summary` file and keep unchanged; the third column is the gene name used in the database file and can be changed to modify the `NAME` field in the HMM profile file; the fourth column and onwards are ignored. If you want to add genes, the intermediate file `HMM.summary` (with the same format as `DB.summary`) located in the `oatkdb` temporary folder gives a full list available genes.

Here is the full list of the options, which can be shown with the command `mancdb -h`.

```
  Program: mancdb (curate OatkDB gene HMM profile database)
  Version: 1.0

  Usage: mancdb [options] gene_list tmpdir output
  Optional:
      -f|--force           Overwrite existing files.
      -s|--sort            Sort gene by names.
         --log      STR    Log file (default stdout).
      -h|--help            Print this help message.
      -v|--version         Print version number.

  Example: mancdb gene_list tmpdir angiosperms_pltd_v20230911
```

## About Gene Names

### tRNA gene
The names of tRNA genes follow the form `trn([A-Z])-([ACGU]{3})`, where the single letter captured by the first matching group is the amino acid code and the three-letter block captured by the second matching group is the anticodon. The program tries to rescue a tRNA gene name violating this naming scheme as long as the name can be matched by the regex `^([t|T][r|R][n|N])([a-z]{0,1})([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]{1})[-|_]?[\(|\]]?([ACGUacgu]{3})?[\)|\]]?$`. If the anticodon matching group is missed, the program will try to retrieve it from the `/note` feature in the Genbank record.

### rRNA gene
The names of rRNA genes follow the form `rrn(\d+(\.\d+)?)`, where the decimal number captured by the first matching group is the product ribosomal subunit. For example, `rrn4.5` and `rrn23` for 4.5S and 23S subunits respectively. If an rRNA gene name violates this naming scheme, the program will try to rescue it by checking the `/product` feature in the Genbank record to retrieve the product subunit. For example, the `MT-RNR1` and `MT-RNR2` genes in the mammal mitochondria genome will be renamed to `rrn12` and `rrn16`, respectively. If you want to change the name back, you can directly edit the HMM profile database. 

### Gene name aliases
A gene could have many names due to different naming conventions used by different annotation tools. For example, the prefix of NADH dehydrogenase genes could be `ND`, `NAD`, `NADH`, `NDH`, etc. OatkDB uses an alias file to capture these different naming schemes for identical genes. An alias file is a two-column (the extra columns will be ignored) text file with each line representing a formatting rule. The first column is the standard name, and the second column is a regex used to capture the matched name aliases.

Here are some examples,

*Example 1.* `COX2  ^[Cc][Oo][Xx]?(2|II)$`

In this example, following this rule, `CO2`, `COII`, `COX2`, `COXII`, `coxII`, etc., will all be captured and converted to `COX2`. 

*Example 2.* `ND(1)  ^[Nn][Aa]?[Dd][Hh]?(\d{1})$`

In this example, a (NADH dehydrogenase) gene name with case insensitive prefixes `ND`, `NAD`, `NADH`, and `NDH` followed by a single digit will be captured and converted to the standard name `ND` followed by the single digit. The `(1)` part in `ND(1)` will be replaced by the first matching group in the formatting regex, i.e., the single digit. There could be more than one matching group, e.g., with the standard name `ND(3)-(1)`, the formatting step will replace `(3)` with the third matching group and `(1)` with the first matching group. In this case, the matching regex (i.e., the second column) must contain at least three matching groups. The program will do this compatibility check and will invoke an error if the check fails.

By default, the `alias.txt` file in this repo will be used by the `oatkdb` script for name alias formatting. The script also accepts a user-defined alias file with `-a|--alias` option. The alias-checking procedure processes the file line by line and stops after the first match, so the user should guarantee no overlaps between matching regexes.

It is noteworthy that the standard names in the alias file are not necessarily the names finally used in the database. The formatting step is only to group together genes with different names. The most commonly used name will be the final choice. For example, in the above example 1, if we see `COX2` 100 times and `COII` 101 times, the name `COII` will be used instead of the standard name `COX2` as specified in the rule, and will be counted 201 times.

