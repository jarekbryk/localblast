## README

Here are some instructions on how to set up BLAST on your local machine (as opposed to run it at NCBI on the web). It is useful if you have lots of sequences to BLAST or if you want to manipulate its output for easier processing downstream.

I included a short presentation on BLAST and the syntax of its run, as well as an example query and a database (see details below) so that you could try it yourself. What I do _not_ cover in these instructions is how to prepare a custom BLAST database from your sequences (if you can't rely on ready-made databases at NCBI), but it is covered in section 7 of the description of BLAST databases from NCBI (included here in full, below).

### Download and install the programme

- [Download BLAST Software and Databases Documentation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs&amp;DOC_TYPE=Download)

- [Standalone BLAST Setup for Unix - BLAST® Help - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK52640/)

- [Standalone BLAST Setup for Windows PC - BLAST® Help - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK52637/)

### Prepare the query

The query for this exercise is a list of 32162 sequences of all unique oligonucleotide probes from [Agilent's _D. melanogaster_ gene expression microarray](http://www.genomics.agilent.com/en/Gene-Expression-Model-Organism-Non-Human-Microarrays/Model-Org-Non-Human-GeneEx-Microarrays/) in a multi-fasta format. The query is included in this repository directly and is called `blast_query.tar.gz`. Uncompressed file is 2.4 MB. Download the query and unzip it into a folder (let's call the folder `blast_pratice` for the purpose of this exercise).

### Prepare the database

The database contains list of all known and predicted transcripts from _D. melanogaster_ in the Ensembl format, as downloaded in March 2017. The database is [available to download from Figshare](https://figshare.com/account/projects/23962/articles/5306095) and the compressed file is approximately 100 MB. Download it and unzip into the `blast_practice` folder - it will unzip into its own folder called `Dmel_transcripts_Ensembl`.

I created those to map Agilent probes against known and predicted transcripts to eliminate probes that do not align well. I used this data in the upcoming publication XXX. If I have some time, I may reduce both the queries and the database for easier sharing and downloading.

### BLAST off!

To run the BLAST (assuming it is installed), enter the `blast_practice` folder and run a megablast of the query against the database with the following command:

`blastn -task megablast -db Dmel_transcripts_Ensembl/Dmel_genes_all.fa -query blast_query.txt -dust no -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue pident stitle" -out outputfile.txt`

The various options are explained in the presentation and in the section **BLAST options** below.

### Useful links

- [BLAST® Command Line Applications User Manual - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

- [BLAST+ features - BLAST® Command Line Applications User Manual - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK279668/#usermanual.Custom_output_formats_for_BLA)

- [Extracting data from BLAST databases with blastdbcmd - BLAST® Command Line Applications User Manual - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK279689/#cookbook.Custom_data_extraction_and_form)

- [Options for the command-line applications. - BLAST® Command Line Applications User Manual - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK279675/)

- [BLAST® Help - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK1762/)

- [Taxonomy browser (Mus musculus)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=10090&amp;lvl=3&amp;lin=f&amp;keep=1&amp;srchmode=1&amp;unlock)

- [BLAST Glossary - BLAST® Help - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK62051/)

### BLAST databases available at NCBI

This is copied directly from the NCBI website and it was up to date as of March 6th, 2017.

<pre>
                         The BLAST Databases
                    Last updated on March 6, 2017

This document describes the BLAST databases available on the NCBI FTP site under 
the /blast/db directory. The direct URL is ftp://ftp.ncbi.nlm.nih.gov/blast/db

1. Quick Start
    * Get all numbered files for a database with the same base name:
      Each of these files represents a subset (volume) of that database,
      and all of them are needed to reconstitute the database.
    * After extraction, there is no need to concatenate the resulting files:
      Call the database with the base name, for nr database files, use "-db nr".
    * For easy download, use the update_blastdb.pl script from the blast+ package.
    * Incremental update is not available.

2. General Introduction

BLAST search pages under the Basic BLAST section of the NCBI BLAST home page
(http://blast.ncbi.nlm.nih.gov/) use a standard set of BLAST databases for 
nucleotide, protein, and translated BLAST searches.  These databases are made 
available as compressed archives of pre-formatted form) and can be donwloaded from
the /db directory of the BLAST ftp site (ftp://ftp.ncbi.nlm.nih.gov/blast/db/). 
The FASTA files reside under the /FASTA subdirectory.

The pre-formatted databases offer the following advantages:
    * Pre-formatting removes the need to run makeblastdb;
    * Species-level taxonomy ids are included for each database entry;
    * Databases are broken into smaller-sized volumes and are therefore easier 
      to download;
    * Sequences in FASTA format can be generated from the pre-formatted databases
      by using the blastdbcmd utility;
    * A convenient script (update_blastdb.pl) is available in the blast+ package 
      to download the pre-formatted databases.

Pre-formatted databases must be downloaded using the update_blastdb.pl script or 
via FTP in binary mode. Documentation for this script can be obtained by running 
the script without any arguments; Perl installation is required.

The compressed files downloaded must be inflated with gzip or other decompress 
tools. The BLAST database files can then be extracted out of the resulting tar 
file using the tar utility on Unix/Linux, or WinZip and StuffIt Expander on 
Windows and Macintosh platforms, respectively.  

Large databases are formatted in multiple one-gigabyte volumes, which are named 
using the basename.##.tar.gz convention. All volumes with the same base name are 
required. An alias file is provided to tie individual volumes together so that 
the database can be called using the base name (without the .nal or .pal 
extension). For example, to call the est database, simply use "-db est" option 
in the command line (without the quotes). 

Additional BLAST databases that are not provided in pre-formatted formats may 
be available in the FASTA subdirectory. For other genomic BLAST databases, 
please check the genomes ftp directory at:
    ftp://ftp.ncbi.nlm.nih.gov/genomes/

3. Contents of the /blast/db/ directory

The pre-formatted BLAST databases are archived in this directory. The names of 
these databases and their contents are listed below.

+-----------------------------+------------------------------------------------+
 File Name                    | Content Description                           
+-----------------------------+------------------------------------------------+
16SMicrobial.tar.gz	          | Bacterial and Archaeal 16S rRNA sequences from 
                                BioProjects 33175 and 33117
FASTA/                        | Subdirectory for FASTA formatted sequences
README                        | README for this subdirectory (this file)
Representative_Genomes.*tar.gz| Representative bacterial/archaeal genomes database
cdd_delta.tar.gz              | Conserved Domain Database sequences for use with 
                                stand alone deltablast
cloud/	                      | Subdirectory of databases for BLAST AMI; see
                                http://1.usa.gov/TJAnEt
env_nr.*tar.gz                | Protein sequences for metagenomes
env_nt.*tar.gz                | Nucleotide sequences for metagenomes
est.tar.gz                    | This file requires est_human.*.tar.gz, 
                                est_mouse.*.tar.gz, and est_others.*.tar.gz files 
                                to function. It contains the est.nal alias so that 
                                searches against est (-db est) will include 
                                est_human, est_mouse and est_others. 
est_human.*.tar.gz            | Human subset of the est database from the est
                                division of GenBank, EMBL and DDBJ.
est_mouse.*.tar.gz            | Mouse subset of the est databasae
est_others.*.tar.gz           | Non-human and non-mouse subset of the est database
gss.*tar.gz                   | Sequences from the GSS division of GenBank, 
                                EMBL, and DDBJ
htgs.*tar.gz                  | Sequences from the HTG division of GenBank, EMBL,
                                and DDBJ
human_genomic.*tar.gz         | Human RefSeq (NC_### ##) chromosome records with 
                                gap adjusted concatenated NT_ contigs
nr.*tar.gz                    | Non-redundant protein sequences from GenPept, 
                                Swissprot, PIR, PDF, PDB, and NCBI RefSeq
nt.*tar.gz                    | Partially non-redundant nucleotide sequences from 
                                all traditional divisions of GenBank, EMBL, and DDBJ 
                                excluding GSS,STS, PAT, EST, HTG, and WGS.
other_genomic.*tar.gz         | RefSeq chromosome records (NC_### ##) for non-human
                                organisms
pataa.*tar.gz                 | Patent protein sequences
patnt.*tar.gz                 | Patent nucleotide sequences. Both patent databases
                                are directly from the USPTO, or from the EPO/JPO 
                                via EMBL/DDBJ
pdbaa.*tar.gz                 | Sequences for the protein structure from the 
                                Protein Data Bank
pdbnt.*tar.gz                 | Sequences for the nucleotide structure from the 
                                Protein Data Bank. They are NOT the protein coding
                                sequences for the corresponding pdbaa entries.
refseq_genomic.*tar.gz        | NCBI genomic reference sequences
refseq_protein.*tar.gz        | NCBI protein reference sequences
refseq_rna.*tar.gz            | NCBI Transcript reference sequences
sts.*tar.gz                   | Sequences from the STS division of GenBank, EMBL,
                                and DDBJ
swissprot.tar.gz              | Swiss-Prot sequence database (last major update)
taxdb.tar.gz                  | Additional taxonomy information for the databases 
                                listed here 
                              | providing common and scientific names
tsa_nt.*tar.gz                | Sequences from the TSA division of GenBank, EMBL,
                                and DDBJ
vector.tar.gz                 | Vector sequences from 2010, see Note 2 in section 4.
+-----------------------------+------------------------------------------------+

4. Contents of the /blast/db/FASTA directory

This directory contains FASTA formatted sequence files. The file names 
and database contents are listed below. These files must be unpacked and 
processed through blastdbcmd before they can be used by the BLAST programs. 

+-----------------------+-----------------------------------------------------+
|File Name              | Content Description                                 |
+-----------------------+-----------------------------------------------------+
alu.a.gz                | translation of alu.n repeats
alu.n.gz                | alu repeat elements (from 2003)
drosoph.aa.gz           | CDS translations from drosophila.nt  
drosoph.nt.gz           | genomic sequences for drosophila (from 2003)
env_nr.gz*              | Protein sequences for metagenomes, taxid 408169
env_nt.gz*              | Nucleotide sequences for metagenomes, taxid 408169
est_human.gz*           | human subset of the est database (see Note 1)
est_mouse.gz*           | mouse subset of the est database
est_others.gz*          | non-human and non-mouse subset of the est database
gss.gz*                 | sequences from the GSS division of GenBank, EMBL, 
                          and DDBJ
htgs.gz*                | sequences from the HTG division of GenBank, EMBL, 
                          and DDBJ 
human_genomic.gz*       | human RefSeq (NC_### ##) chromosome records
                          with gap adjusted concatenated NT_ contigs 
igSeqNt.gz              | human and mouse immunoglobulin variable region 
                          nucleotide sequences
igSeqProt.gz            | human and mouse immunoglobulin variable region 
                          protein sequences
mito.aa.gz              | CDS translations of complete mitochondrial genomes
mito.nt.gz              | complete mitochondrial genomes
nr.gz*                  | non-redundant protein sequence database with entries
                           from GenPept, Swissprot, PIR, PDF, PDB, and RefSeq
nt.gz*                  | nucleotide sequence database, with entries from all 
                          traditional divisions of GenBank, EMBL, and DDBJ; 
                          excluding bulk divisions (gss, sts, pat, est, htg) 
                          and wgs entries. Partially non-redundant.
other_genomic.gz*       | RefSeq chromosome records (NC_### ##) for organisms 
                          other than human
pataa.gz*               | patent protein sequences
patnt.gz*               | patent nucleotide sequences. Both patent sequence 
                          files are from the USPTO, or EPO/JPO via EMBL/DDBJ
pdbaa.gz*               | protein sequences from pdb protein structures
pdbnt.gz*               | nucleotide sequences from pdb nucleic acid 
                          structures. They are NOT the protein coding 
                          sequences for the corresponding pdbaa entries.
sts.gz*                 | database for sequence tag site entries 
swissprot.gz*           | swiss-prot database (last major release)
vector.gz               | vector sequences from 2010. (See Note 2)
yeast.aa.gz             | protein translations from yeast.nt
yeast.nt.gz             | yeast genomes (from 2003)
+-----------------------+---------------------------------------------------+
NOTE: 
(1) NCBI does not provide the complete est database in FASTA format. One 
    needs to get all three subsets (est_human, est_mouse, and est_others
    and concatenate them into the complete est fasta database).
(2) For screening for vector contamination, use the UniVec database:
    ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/
 *  marked files have pre-formatted counterparts. 

5. Database updates

The BLAST databases are updated regularly. There is no established incremental
update scheme. We recommend downloading the complete databases regularly to 
keep their content current.

6. Non-redundant defline syntax

The non-redundant databases are nr, nt and pataa. Identical sequences are 
merged into one entry in these databases. To be merged two sequences must
have identical lengths and every residue at every position must be the 
same.  The FASTA deflines for the different entries that belong to one 
record are separated by control-A characters invisible to most 
programs. In the example below both entries Q57293.1 and AAB05030.1
have the same sequence, in every respect:

>Q57293.1 RecName: Full=Fe(3+) ions import ATP-binding protein FbpC ^AAAB05030.1 afuC 
[Actinobacillus pleuropneumoniae] ^AAAB17216.1 afuC [Actinobacillus pleuropneumoniae]
MNNDFLVLKNITKSFGKATVIDNLDLVIKRGTMVTLLGPSGCGKTTVLRLVAGLENPTSGQIFIDGEDVTKSSIQNRDIC
IVFQSYALFPHMSIGDNVGYGLRMQGVSNEERKQRVKEALELVDLAGFADRFVDQISGGQQQRVALARALVLKPKVLILD
EPLSNLDANLRRSMREKIRELQQRLGITSLYVTHDQTEAFAVSDEVIVMNKGTIMQKARQKIFIYDRILYSLRNFMGEST
ICDGNLNQGTVSIGDYRFPLHNAADFSVADGACLVGVRPEAIRLTATGETSQRCQIKSAVYMGNHWEIVANWNGKDVLIN
ANPDQFDPDATKAFIHFTEQGIFLLNKE

Individual sequences are now identifed simply by their accession.version.  

For databases whose entries are not from official NCBI sequence databases, 
such as Trace database, the gnl| convention is used. For custom databases, 
this convention should be followed and the id for each sequence must be 
unique, if one would like to take the advantage of indexed database, 
which enables specific sequence retrieval using blastdbcmd program included 
in the blast executable package.  One should refer to documents 
distributed in the standalone BLAST package for more details. 

7. Formatting a FASTA file into a BLASTable database

FASTA files need to be formatted with makeblastdb before they can be used in local 
blast search. For those from NCBI, the following makeblastdb commands are
recommended:

For nucleotide fasta file:   makeblastdb -in input_db -dbtype nucl -parse_seqids
For protein fasta file:      makeblastdb -in input_db -dbtype prot -parse_seqids

In general, if the database is available as BLAST database, it is better to use the 
preformatted database.

8. Technical Support

Questions and comments on this document and NCBI BLAST related questions 
should be sent to the blast-help group at:
      blast-help@ncbi.nlm.nih.gov

For information about other NCBI resources/services, please send email to 
NCBI User Service at:
      info@ncbi.nlm.nih.gov
</pre>

### BLAST options

**DESCRIPTION**

Nucleotide-Nucleotide BLAST 2.2.31+

**USAGE**

`blastn [-h] [-help] [-import_search_strategy filename] [-export_search_strategy filename] [-task task_name] [-db database_name] [-dbsize num_letters] [-gilist filename] [-seqidlist filename] [-negative_gilist filename] [-entrez_query entrez_query] [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm] [-subject subject_input_file] [-subject_loc range] [-query input_file] [-out output_file] [-evalue evalue] [-word_size int_value] [-gapopen open_penalty] [-gapextend extend_penalty] [-perc_identity float_value] [-qcov_hsp_perc float_value] [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value] [-xdrop_gap_final float_value] [-searchsp int_value] [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy] [-min_raw_gapped_score int_value] [-template_type type] [-template_length int_value] [-dust DUST_options] [-filtering_db filtering_database] [-window_masker_taxid window_masker_taxid] [-window_masker_db window_masker_db] [-soft_masking soft_masking] [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value] [-best_hit_score_edge float_value] [-window_size int_value] [-off_diagonal_range int_value] [-use_index boolean] [-index_name string] [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format] [-show_gis] [-num_descriptions int_value] [-num_alignments int_value] [-line_length line_length] [-html] [-max_target_seqs num_sequences] [-num_threads int_value] [-remote] [-version]`

**OPTIONAL ARGUMENTS**

- `-h` Print USAGE and DESCRIPTION; ignore all other parameters
- `-help` Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
- `-version` Print version number; ignore other arguments

### Input query options

- `-query <File_In>` Input file name; Default = '-'
- `-query_loc <String>` Location on the query sequence in 1-based offsets (Format: start-stop)
- `-strand <String, 'both', 'minus', 'plus'>` Query strand(s) to search against database/subject; Default = 'both'

### General search options
 
- `-task <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn' >` Task to execute; Default = `megablast`
- `-db <String>` BLAST database name; Incompatible with: `subject, subject_loc`
- `-out <File_Out>` Output file name; Default = - `-`
- `-evalue <Real>` Expectation value (E) threshold for saving hits; Default = `10`
- `-word_size <Integer, >=4>` Word size for wordfinder algorithm (length of best perfect match)
- `-gapopen <Integer>` Cost to open a gap
- `-gapextend <Integer>` Cost to extend a gap
- `-penalty <Integer, <=0>` Penalty for a nucleotide mismatch
- `-reward <Integer, >=0>` Reward for a nucleotide match
- `-use_index <Boolean>` Use MegaBLAST database index; Default = `false`
- `-index_name <String>` MegaBLAST database index name

### BLAST-2-Sequences options
 
- `-subject <File_In>` Subject sequence(s) to search; Incompatible with: `db, gilist, seqidlist, negative_gilist, db_soft_mask, db_hard_mask`
- `-subject_loc <String>` Location on the subject sequence in 1-based offsets (Format: start-stop); Incompatible with: `db, gilist, seqidlist, negative_gilist, db_soft_mask, db_hard_mask, remote`

### Formatting options

- `-outfmt <String>` alignment view options:
	- `0` = pairwise,
	- `1` = query-anchored showing identities,
	- `2` = query-anchored no identities,
	- `3` = flat query-anchored, show identities,
	- `4` = flat query-anchored, no identities,
	- `5` = XML Blast output,
	- `6` = tabular,
	- `7` = tabular with comment lines,
	- `8` = Text ASN.1,
	- `9` = Binary ASN.1,
	- `10` = Comma-separated values,
	- `11` = BLAST archive format (ASN.1),
	- `12` = JSON Seqalign output,
	- `13` = JSON Blast output,
	- `14` = XML2 Blast output
 
Options `6`, `7`, and `10` can be additionally configured to produce a custom format specified by space delimited format specifiers. The supported format specifiers are:

- `-outfmt <String>`
	- `qseqid` means Query Seq-id
	- `qgi` means Query GI
	- `qacc` means Query accesion
	- `qaccver` means Query accesion.version
	- `qlen` means Query sequence length
	- `sseqid` means Subject Seq-id
	- `sallseqid` means All subject Seq-id(s), separated by a ';'
	- `sgi` means Subject GI
	- `sallgi` means All subject GIs
	- `sacc` means Subject accession
	- `saccver` means Subject accession.version
	- `sallacc` means All subject accessions
	- `slen` means Subject sequence length
	- `qstart` means Start of alignment in query
	- `qend` means End of alignment in query
	- `sstart` means Start of alignment in subject
	- `send` means End of alignment in subject
	- `qseq` means Aligned part of query sequence
	- `sseq` means Aligned part of subject sequence
	- `value` means Expect value
	- `bitscore` means Bit score
	- `score` means Raw score
	- `length` means Alignment length
	- `pident` means Percentage of identical matches
	- `nident` means Number of identical matches
	- `mismatch` means Number of mismatches
	- `positive` means Number of positive-scoring matches
	- `gapopen` means Number of gap openings
	- `gaps` means Total number of gaps
	- `ppos` means Percentage of positive-scoring matches
	- `frames` means Query and subject frames separated by a '/'
	- `qframe` means Query frame
	- `sframe` means Subject frame
	- `btop` means Blast traceback operations (BTOP)
	- `staxids` means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
	- `sscinames` means unique Subject Scientific Name(s), separated by a ';'
	- `scomnames` means unique Subject Common Name(s), separated by a ';'
	- `blastnames` means unique Subject Blast Name(s), separated by a ';' (in alphabetical order)
	- `sskingdoms` means unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)
	- `stitle` means Subject Title
	- `salltitles` means All Subject Title(s), separated by a '<>'
	- `sstrand` means Subject Strand
	- `qcovs` means Query Coverage Per Subject
	- `qcovhsp` means Query Coverage Per HSP
	- `qcovus` is a measure of Query Coverage that counts a position in a subject sequence for this measure only once. The second time the position is aligned to the query is not counted towards this measure.

When - `-outfmt` options are not specified, the default value is: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`, which is equivalent to the keyword `std`.

- `-show_gis` Show NCBI GIs in deflines
- `-num_descriptions <Integer, >=0>` Number of database sequences to show one-line descriptions for; Not applicable for `outfmt >4`; Default = `500`; Incompatible with: `max_target_seqs`
- `-num_alignments <Integer, >=0>` Number of database sequences to show alignments for; Default = `250`; Incompatible with: `max_target_seqs`
- `-line_length <Integer, >=1>` Line length for formatting alignments; Not applicable for `outfmt >4`; Default = `60`
- `-html` Produce HTML output

### Query filtering options

- `-dust <String>` Filter query sequence with DUST (Format: 'yes', 'level window linker', or 'no' to disable); Default = `20 64 1`
- `-filtering_db <String>` BLAST database containing filtering elements (i.e.: repeats)
- `-window_masker_taxid <Integer>` Enable WindowMasker filtering using a Taxonomic ID
- `-window_masker_db <String>` Enable WindowMasker filtering using this repeats database.
- `-soft_masking <Boolean>` Apply filtering locations as soft masks
 Default = `true'
- `-lcase_masking
 Use lower case filtering in query and subject sequence(s)?

### Restrict search or results

- `-gilist <String>` Restrict search of database to list of GI's; Incompatible with: `negative_gilist, seqidlist, remote, subject, subject_loc`
- `-seqidlist <String>` Restrict search of database to list of SeqId's; Incompatible with: `gilist, negative_gilist, remote, subject, subject_loc`
- `-negative_gilist <String>` Restrict search of database to everything except the listed GIs; Incompatible with: `gilist, seqidlist, remote, subject, subject_loc`
- `-entrez_query <String>` Restrict search with the given Entrez query; Requires: `remote`
- `-db_soft_mask <String>` Filtering algorithm ID to apply to the BLAST database as soft masking; Incompatible with: `db_hard_mask, subject, subject_loc`
- `-db_hard_mask <String>` Filtering algorithm ID to apply to the BLAST database as hard masking; Incompatible with: `db_soft_mask, subject, subject_loc`
- `-perc_identity <Real, 0..100>` Percent identity
- `-qcov_hsp_perc <Real, 0..100>` Percent query coverage per hsp
- `-max_hsps <Integer, >=1>` Set maximum number of HSPs per subject sequence to save for each query
- `-culling_limit <Integer, >=0>` If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit; Incompatible with: `best_hit_overhang, best_hit_score_edge`
- `-best_hit_overhang <Real, (>0 and <0.5)>` Best Hit algorithm overhang value (recommended value: 0.1); Incompatible with: `culling_limit`
- `-best_hit_score_edge <Real, (>0 and <0.5)>` Best Hit algorithm score edge value (recommended value: 0.1); Incompatible with: `culling_limit`
- `-max_target_seqs <Integer, >=1>` Maximum number of aligned sequences to keep; Not applicable for `outfmt <=4`; Default = `500`; Incompatible with: `num_descriptions, num_alignments`

### Discontiguous MegaBLAST options

- `-template_type <String, coding, coding_and_optimal, optimal>` Discontiguous MegaBLAST template type; Requires: `template_length`
- `-template_length <Integer, Permissible values: 16, 18, 21>` Discontiguous MegaBLAST template length; Requires: `template_type`

### Statistical options

- `-dbsize <Int8>` Effective length of the database 
- `-searchsp <Int8, >`=0>` Effective length of the search space
- `-sum_stats <Boolean>` Use sum statistics

### Search strategy options

- `-import_search_strategy <File_In>` Search strategy to use; Incompatible with:  `export_search_strategy`
- `-export_search_strategy <File_Out>` File name to record the search strategy used; Incompatible with: `import_search_strategy`

### Extension options
 
- `-xdrop_ungap <Real>` X-dropoff value (in bits) for ungapped extensions
- `-xdrop_gap <Real>` X-dropoff value (in bits) for preliminary gapped extensions
- `-xdrop_gap_final <Real>` X-dropoff value (in bits) for final gapped alignment
- `-no_greedy` Use non-greedy dynamic programming extension
- `-min_raw_gapped_score <Integer>` Minimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages
- `-ungapped` Perform ungapped alignment only
- `-window_size <Integer, >=0>` Multiple hits window size, use 0 to specify 1-hit algorithm
- `-off_diagonal_range <Integer, >=0>` Number of off-diagonals to search for the 2nd hit, use `0` to turn off; Default = `0`

### Miscellaneous options

- `-parse_deflines` Should the query and subject defline(s) be parsed?
- `-num_threads <Integer, >=1>` Number of threads (CPUs) to use in the BLAST search; Default = `1`; Incompatible with: `remote`
- `-remote` Execute search remotely?; * Incompatible with: `gilist, seqidlist, negative_gilist, subject_loc, num_threads`