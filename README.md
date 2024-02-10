# FastQHandler
A collection of Python scripts designed for the efficient management of various FASTQ file formats.

![FastQHandler Logo](https://github.com/OZTaekOppa/FastQHandler/blob/main/images/FastQHandler_Logo.png)


## Brief Background
**FastQHandler**, created by Hyungtaek Jung and the team at the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au/) at [The Australian National University](https://www.anu.edu.au/), is a Python script suite for efficient FASTQ file management. It boasts 17 work modules to ease input/output processes, covering various aspects of FASTQ data analysis including post-processing and format conversion. Optimized for life science datasets, **FastQHandler** is a CLI application tested across different FASTQ formats. Users should note that processing large datasets may require substantial computational resources on Linux, HPC or Cloud platforms.


## Citation
Hyungtaek Jung, Kirat Alreja, Kosar Hooshmand, Hadi Nazem-Bokaee, Hardip Patel: **FastQHandler**: An easy Python-based toolset for handling fasta files, [PLoS Comp Biol Submitted](https://www.biorxiv.org/XXXX).


## Contents:
+ [STABLE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#stable-version-101)
+ [INSTALLATION](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#installation)
+ [LICENSE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#license) 
+ [GETTING STARTED](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#getting-started)
+ [FAQ](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#faq)
+ [WIKI PAGE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#wiki-page)
+ [AUTHORS](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#authors)
+ [COPYRIGHT](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#copyright)


## STABLE (version 1.0.1)
- Release date: February 2024
- **FastQHandler** is a standalone Python application equipped with 17 sub-modules for interactive FASTQ file manipulation, available as open-source (see [LICENSE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#license)).


## INSTALLATION
- Access the program via [GitHub](https://github.com/OZTaekOppa/FastQHandler)
- Installation options include Bioconda or Python pip. For support, refer to [Issues on GitHub](https://github.com/OZTaekOppa/FastQHandler/issues).
- **FastQHandler** requires no separate installation process.
- Just clone this repository, and run
```
git clone https://github.com/OZTaekOppa/FastQHandler/
python3 {path}/fastqhandler.py
```

## LICENSE
**FastQHandler** is available under the MIT license and incorporates various open-source software. For detailed information on the integrated Python packages, modules, and libraries, and their specific applications within **FastQHandler**, please refer to the [manuscript](https://www.biorxiv.org/XXXX)


## GETTING STARTED
**FastQHandler** is developed primarily in Python 3.9+ and Biopython and features 17 modules. It facilitates data input and output through a Command-Line Interface (CLI), ensuring smooth end-to-end file handling. To optimize the use of **FastQHandler**, users should prepare all necessary input files, such as FASTQ and TXT formats, in advance.

![FASTAhandler Workflow](https://github.com/OZTaekOppa/FastQHandler/blob/main/images/FastaHandler_Workflow.png)

### Tested FASTQ Files
- Illumina Paired-End: EBI [SRRSRR26400271 NovoSeq 6000](https://www.ebi.ac.uk/ena/browser/view/SRX22106079) and [SRR26087471 Genome Analyzer](https://www.ebi.ac.uk/ena/browser/view/SRR26087471)
- PacBio: EBI [SRR10326407 Sequel](https://www.ebi.ac.uk/ena/browser/view/SRR10326407)
- Oxford Nanopore: EBI [SRR12145915 MinION](https://www.ebi.ac.uk/ena/browser/view/SRR12145915)

### General Usage
```
FastQHandler: Fastq File Manipulation Toolkit
version 1.0.1

Usage: python3 fastqhandler.py <module> <parameters>

Modules:
Multi2Single    		| m2s   	Convert a multi-fasta (multiline) into a single-line fasta.
RenameId        		| rid   	Rename prefix IDs and headers.
PrefixRename    		| prn   	Rename prefix IDs and headers with a user’s input.
PrefixSelectRename      	| psr   	Rename prefix IDs and headers with a user’s input (Only).
IdExtract       		| idx   	Extract matched IDs and their corresponding sequences.
IdExtractLocation       	| iel   	Extract matched IDs, locations and their corresponding sequences.
IdExtractLocationMultiple       | iem   	Extract matched IDs, locations and their corresponding sequences (Multiple).
ReverseComplement       	| rcp   	Make a reverse complement sequence.
FindCountDuplication    	| fcd   	Find and count the duplicated IDs and sequences.
RemoveDuplication       	| rvd   	Remove the duplicated IDs and sequences.
SubsetFasta     		| ssf   	Make a subset of data with a sequence length filter.
ExtractPattern  		| xpt   	Make a subset of data with find, filter and extract.
EachFastaStats  		| efs   	Generate each line fasta statistic for a multi-line fasta.
AllFastaStats   		| afs   	Generate a summary of multi-line fasta statistics.
MultipleFastaStats      	| mfs   	Generate a summary of multi-line fasta statistics (Multiple).
ConcatenateFasta        	| ccf   	Make a concatenated fasta file (Multiple).
TranslateSequence       	| tls   	Find the translated sequences as a protein and open reading frames (ORFs).

Use <module> --help for module usage.
```

### fqread_stats (qrs)
- To generate fastq statistics for a fastq file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A summary of fastq with its diverse statistics (e.g. sequence length, max length, N50, and more).

Example usage
```
python3 fqread_stats.py --input-seq test_dna.fq --out output_stats.csv --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fqread_stats.py:  Call fqread_stats module
	1. python3 fqread_stats.py --help: Check help menu
		+ --input-seq: Indicate an input multi-line fastq file and its path
		+ --out: Indicate an output csv file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqread_stats_unlim (qsu)
- To generate fastq statistics for multiple fastq files
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multiple fastq files.
	+ Output: A summary of fastq with its diverse statistics (e.g. sequence length, max length, N50, and more).

Example usage
```
python3 fqread_stats_unlim.py --input-seqs test_dna1.fq test_dna2.fq test_dna3.fq test_dna4.fq test_dna5.fq --out ./ --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fqread_stats_unlim.py:  Call fqread_stats_unlim module
	1. python3 fqread_stats_unlim.py --help: Check help menu
		+ --input-seqs: Indicate multiple fastq files and theirs paths
		+ --out: Indicate an output folder to save a csv file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fq2fa_convert (qac)
- To convert a fastq file as a fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A converted fasta file with the same fastq header.

Example usage
```
python3 fq2fa_convert.py --input-seq test_dna.fq --out output_convert.fq --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fq2fa_convert.py:  Call fq2fa_convert module
	1. python3 fq2fa_convert.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate a converted output fasta file with the same fastq header
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqlen_filter (qlf) and fqlen_filter_only (qlo)
- To extract filtered IDs and their corresponding sequences (focused on sequence length)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences based on a sequence length filter.

Example usage
```
python3 fqlen_filter.py --input-seq test_dna.fq --filter-len 100 --out ./ --t 1 --mem 2
```

+ The script will save sequences that pass and fail the filtering criteria into filter_passed.fq and filter_failed.fq, respectively.
+ To generate a single outcome of filter_passed.fq, please use **fqlen_filter_only (qlo)** script with a specific output file name.

- Parameter explanation
	1. python 3: Call python 3
	1. fqlen_filter.py:  Call fqlen_filter module
	1. python3 fqlen_filter.py --help: Check help menu
	--input-seq: Indicate an input fastq file and its path
	--filter-len: Indicate a sequence length to filter out
	--out: Indicate an output ID filtered and extracted single-line fastq file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)



### fqqual_filter (qqf) and fqqual_filter_only (qqo)
- To extract filtered IDs and their corresponding sequences (focused on a read quality)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences based on a read quality.

Example usage
```
python3 fqqual_filter.py --input-seq test_dna.fq --filter-qual 15 --out ./ --t 1 --mem 2
```

+ The script will save sequences that pass and fail the filtering criteria into filter_passed.fq and filter_failed.fq, respectively.
+ To generate a single outcome of filter_passed.fq, please use **fqqual_filter_only (qqo)** script with a specific output file name.

  
- Parameter explanation
	1. python 3: Call python 3
	1. fqqual_filter.py:  Call fqqual_filter module
	1. python3 fqqual_filter.py --help: Check help menu
	--input-seq: Indicate an input fastq file and its path
	--filter-qual: Indicate a read quality to filter out
	--out: Indicate an output ID filtered and extracted single-line fastq file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)


### fqlenqual_filter (q2f) and fqlenqual_filter_only (q2o)
- To extract filtered IDs and their corresponding sequences (focused on a sequence length and read quality)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences based on a sequence length and read quality.

Example usage
```
python3 fqlenqual_filter.py --input-seq test_dna.fq --filter-qual 15 filter-len 100 --out ./ --t 1 --mem 2
```

+ The script will save sequences that pass and fail the filtering criteria into filter_passed.fq and filter_failed.fq, respectively.
+ To generate a single outcome of filter_passed.fq, please use **fqlenqual_filter_only (q2o)** script with a specific output file name.

- Parameter explanation
	1. python 3: Call python 3
	1. fqlenqual_filter.py:  Call fqlenqual_filter module
	1. python3 fqlenqual_filter.py --help: Check help menu
	--input-seq: Indicate an input fastq file and its path
	--filter-qual: Indicate a read quality to filter out
	--filter-len: Indicate a sequence length to filter out
	--out: Indicate an output ID filtered and extracted single-line fastq file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)


### fqtrim (qtr)
- To extract trimmed IDs and their corresponding sequences (focused on a sequence length; right, left or both)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences based on a sequence length.

Example usage
```
python3 fqtrim.py --input-seq test_dna.fq --trim-g 5 --trim-a 10 --trim-b 10 --out output_trimmed.fq --t 1 --mem 2
```

+ Please indicate one of the parameters to trim (-g, -a or -b). If not, a default parameter, "-b 10", will be applied.

- Parameter explanation
	1. python 3: Call python 3
	1. fqtrim.py:  Call fqtrim module
	1. python3 fqtrim.py --help: Check help menu
	--input-seq: Indicate an input fastq file and its path
	--trim-g: Indicate a read length to trim from the left (5'end)
	--trim-a: Indicate a read length to trim from the right (3'end)
	--trim-b: Indicate a read length to trim from both (5' and 3'end)
	--out: Indicate an output ID trimmed and extracted single-line fastq file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)


### fqseq_find (qsf)
- To extract matched IDs and their locations (focused on a matched sequence)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file and a sequence to find.
	+ Output: A csv with extracted IDs and their corresponding locations (start and end).

Example usage
```
python3 fqseq_find.py --input-seq test_dna.fq --find-seq ATGCCGGTGAC --out output_found.csv --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. fqseq_find.py:  Call fqseq_find module
	1. python3 fqseq_find.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --find-seq: Indicate a sequence to find a pattern
  		+ --out: Indicate an output ID matched and extracted fastq file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqseq_repeat_find (qrf)
- To extract matched IDs and their counts (focused on a matched and repeated sequence)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file and a sequence to find.
	+ Output: A matched fastq and a csv with extracted IDs and their counts.
	+ Example file: [rpt_seq.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/rpt_seq.txt) in the "example_data" folder.

Example usage
```
python3 fqseq_repeat_find.py --input-seq test_dna.fq --find-rpt-seq rpt_seq.txt --out ./ --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. fqseq_repeat_find.py:  Call fqseq_repeat_find module
	1. python3 fqseq_repeat_find.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --find-rpt-seq: Indicate a repeat sequence file and its path to find
  		+ --out: Indicate an output ID matched, count and extracted fastq file and its path (both cvs and fastq)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqid_extract (qix)
- To extract matched IDs and their corresponding sequences.
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file and a list of IDs.
	+ Output: A fastq with extracted IDs and their corresponding sequences.
 	+ Example file: [id_extract.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/id_extract.txt) in the "example_data" folder.

Example usage
```
python3 fqid_extract.py --input-seq test_dna.fq --find-id id_extract.txt --out output_extracted.fasta --t 1 --mem 2
```
	
- Parameter explanation
	1. python 3: Call python 3
	1. fqid_extract.py:  Call fqid_extract module
	1. python3 fqid_extract.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --find-id: Indicate an input ID and header (without "@") text file and its path
		+ --out: Indicate an output ID matched and extracted fastq file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqid_rename (qin)
- To rename prefix IDs and headers from a fastq file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A fastq with a new prefix ID name.


Example usage
```
python3 fqid_rename.py --input-seq test_dna.fq --new-id FunNGS --out output_renamed.fq --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. fqid_rename.py:  Call fqid_rename module
	1. python3 fqid_rename.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --new-id: Indicate a new prefix ID/header name (accept both integer and strings but no space)
		+ --out: Indicate an output renamed fastq file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqsubset (qsb)
- To make a subset of data counting from beginning, end or randomly
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A subsetted fastq.

Example usage
```
python3 fqsubset.py --input-seq test_dna.fq --sbs-g 100 --sbs-a 100 --sbs-r 100 --out output_subset.fq --t 1 --mem 2
```

+ Please indicate one of the parameters to subset (-g, -a or -r). If not, a default parameter, "-r 100", will be applied.

- Parameter explanation
	1. python 3: Call python 3
	1. fqsubset.py:  Call fqsubset module
	1. python3 fqsubset.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
	  	+ --trim-g: Indicate a read count to subset from the beginning 
		+ --trim-a: Indicate a read count to subset from the end 
		+ --trim-r: Indicate a read count randomly
		+ --out: Indicate an output fastq file after subsetting
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqsubset_pe (qsp)
- To make a subset of paired-end data (e.g. Illumina) counting from beginning, end or randomly
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A paired-end fastq file.
	+ Output: A subsetted fastq.

Example usage
```
python3 fqsubset_pe.py --input-seq1 test_r1.fq --input-seq2 test_r2.fq --sbs-g 100 --sbs-a 100 --sbs-r 100 --out ./ --t 1 --mem 2
```

+ Please indicate one of the parameters to subset (-g, -a or -r). If not, a default parameter, "-r 100", will be applied.
+ The script will save sequences that two read pairs the subsetting criteria into subset_pe1.fq and subset_pe2.fq, respectively.

- Parameter explanation
	1. python 3: Call python 3
	1. fqsubset.py:  Call fqsubset module
	1. python3 fqsubset.py --help: Check help menu
		+ --input-seq1: Indicate an input r1.fq file and its path
  		+ --input-seq2: Indicate an input r2.fq file and its path
	  	+ --trim-g: Indicate a read count to subset from the beginning 
		+ --trim-a: Indicate a read count to subset from the end 
		+ --trim-r: Indicate a read count randomly
		+ --out: Indicate an output folder for read pairs of fastq after subsetting
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqid_remove (qir)
- To remove the same prefix IDs and headers from a fastq file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A fastq with a summary csv.


Example usage
```
python3 fqid_remove.py --input-seq test_dna.fq --out ./ --t 1 --mem 2
```

+ The script will remove the same prefix IDs and headers if they are duplicated.
  
- Parameter explanation
	1. python 3: Call python 3
	1. fqid_remove.py:  Call fqid_remove module
	1. python3 fqid_remove.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate an output of non-duplicated fastq, csv and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqseq_remove (qsr)
- To remove the same sequence despite the discrepancy of IDs and headers from a fastq file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A fastq with a summary csv.


Example usage
```
python3 fqseq_remove.py --input-seq test_dna.fq --out ./ --t 1 --mem 2
```

+ If they are duplicated, the script will remove the same sequence despite the discrepancy of IDs and headers.
  
- Parameter explanation
	1. python 3: Call python 3
	1. fqseq_remove.py:  Call fqseq_remove module
	1. python3 fqseq_remove.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate an output of non-duplicated fastq, csv and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqdup_count (qdc)
- To find and count the same IDs/headers and the same sequence despite the discrepancy of IDs/headers from a fastq file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A fastq with a summary csv.
		- dupid_count.csv: The same IDs and headers
  		- dupsep_count.csv: The same sequence despite the discrepancy of IDs and headers


Example usage
```
python3 fqdup_count.py --input-seq test_dna.fq --out ./ --t 1 --mem 2
```

+ If they are duplicated (IDs, headers and sequences), the script will find and count the IDs/headers and their occurrence.
  
- Parameter explanation
	1. python 3: Call python 3
	1. fqdup_count.py:  Call fqdup_count module
	1. python3 fqdup_count.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate an output of csv and its path 
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)



### fqrev_complt (rcp)
- To make a reverse complement sequence
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A reverse complement fastq.

Example usage
```
python3 fqrev_complt.py --input-seq test_dna.fq --out output_revcomplt.fq --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. fqrev_complt.py:  Call fqrev_complt module
	1. python3 fqrev_complt.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate a reverse complement converted output fastq file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)




### multi2single (m2s)
- To convert a multi-fasta (multiline) into a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A multiline fasta file.
	+ Output: A single-line fasta.
	+ Example file: [mltseq2sl](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/mltseq2sl.fa) in the "example_data" folder. 

Example usage
```
python3 multi2single.py --input-seq test_dna.fasta --out test_output_sl.fasta --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. multi2singleline.py:  Call multi2single module
	1. python3 multi2singleline.py --help: Check help menu
		+ --input-seq: Indicate an input multi-fasta (multiline) fasta file and its path
		+ --out: Indicate an output single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only only with Gb size)



### prfxrename (prn)
- To rename prefix IDs and headers from a single-line fasta with a user's input text file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a tap-separated text file (e.g. old_IDs	new_IDs)
	+ Output: A single-line fasta with a new prefix ID name based on the user's input text file. Unmatched IDs will be also produced with its original IDs.
	+ Example file: [header_id_only.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/header_id_only.txt) in the "example_data" folder.


Example usage
```
python3 prfxrename.py --input-seq test_dna.fasta --input-id new_ids.txt --out output_reN.fasta --t 1 --mem 2
```
+ Use the multi2single module first if your input FASTA file isn't in single-line format.
+ The script is ideal for renaming specific parts of IDs/headers in FASTA files, such as those from genome assemblies acquired from public databases like NCBI and EBI, based on user input.

- Parameter explanation
	1. python 3: Call python 3
	1. prfxrename.py:  Call prfxrename module
	1. python3 renameid.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-id: Indicate a new tap-separated prefix ID/header name file (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### prfxselrename (psr)
- To rename prefix IDs and headers from a single-line fasta with a user's input text file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a tap-separated text file (e.g. old_IDs	new_IDs)
	+ Output: A single-line fasta with a new prefix ID name based on the user's input text file. Unmatched IDs will be discarded.
 	+ Example file: [header_id_only.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/header_id_only.txt) in the "example_data" folder.


Example usage
```
python3 prfxselrename.py --input-seq test_dna.fasta --input-id new_ids.txt --out output_reN.fasta --t 1 --mem 2
```
+ Use the multi2single module first if your input FASTA file isn't in single-line format.
+ The script is effective for selectively renaming IDs/headers in FASTA files (with a user's specific input), such as genome assemblies from databases, where you want to exclude and not rename certain elements like scaffolds.

- Parameter explanation
	1. python 3: Call python 3
	1. prfxselrename.py:  Call prfxselrename module
	1. python3 renameid.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-id: Indicate a new tap-separated prefix ID/header name file (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)



### idextloct (iel)
- To extract matched IDs, locations and their corresponding sequences (focused on a single ID)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.
	+ Example file: [header_id.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/header_id.txt) and [idextloct.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/idextloct.fa) in the "example_data" folder.

Example usage
```
python3 idextloct.py --input-seq test_dna.fasta --header-id test3_3%week --start 2 --end 10 --out output_test.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. idextractlocation.py:  Call idextloct module
	1. python3 idextloct.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --header-id: Indicate an input ID and header (without ">") name and pattern
		+ --start: Indicate a start position to extract (please use -1 value due to the python index)
		+ --end: Indicate en end position to extract (please use -1 value due to the python index)
		+ --out: Indicate an output ID matched and extracted single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)
	
### idextloctmlt (iem)
- To extract matched IDs, locations and their corresponding sequences (focused on multiple IDs)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.
 	+ Example file: [hdrmulti_id.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/hdrmulti_id.txt), [idextract.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/idextract.fa) and [idextloct.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/idextloct.fa) in the "example_data" folder.

Example usage
```
python3 idextloctmlt.py --input-seq test_dna.fasta --input-extract input_extract.txt --out output_extest.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. idextractlocamulti.py:  Call idextloctmlt module
	1. python3 idextractlocamulti.py --help: Check help menu
	--input-seq: Indicate an input single-line fasta file and its path
	--input-extract: Indicate an input ID and header (without ">") text file (a tap-separated) and its path including start and end positions (please use -1 value due to the python index)
	--out: Indicate an output ID matched and extracted single-line fasta file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)



### findcntdupl (fcd)
- To find and count the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A text with a duplication number.
	+ Example file: [dupids.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/dupids.fa) in the "example_data" folder.

Example usage
```
python3 findcntdupl.py --input-seq test_dna2.fasta --out output_files.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. findcountdupl.py:  Call findcntdupl module
	1. python3 findcntdupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --out: Indicate an output text file after finding and counting the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### removedupl (rvp)
- To remove the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A single-line fasta after removing updated IDs and sequences.
 	+ Example file: [rmvdup.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/rmvdup.fa) in the "example_data" folder.

Example usage
```
python3 removedupl.py --input-seq test_dna2.fasta --outfasta output_testdupl.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. removedupl.py:  Call removedupl module
	1. python3 removedupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --outfasta: Indicate an output fasta file after removing the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### subsetfa (ssf)
- To make a subset of data with a sequence length filter option
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A subseted single-line fasta.
 	+ Example file: [subsetfa.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/subsetfa.fa) in the "example_data" folder.

Example usage
```
python3 subsetfa.py --input-seq test_mRNA1.fasta --filter 50 --out output_subset.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. subsetfasta.py:  Call subsetfa module
	1. python3 subsetfa.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --filter: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (integer only)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### extractptrn (xpt)
- To make a subset of data with find, filter and extract options
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of patterns.
	+ Output: A extracted single-line fasta with IDs and their corresponding sequences.
 	+ Example file: [find_ptrn.txt](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/find_ptrn.txt), [find_ptrn_mlt.txt
](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/find_ptrn_mlt.txt), [find_ptrn_mlt_seq.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/find_ptrn_mlt_seq.fa), and [find_ptrn_seq.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/find_ptrn_seq.fa) in the "example_data" folder.

Example usage
```
python3 extractptrn.py --input-seq test_dna1.fasta --input-pattern seq_pattern.txt --input-length 45 --out output_pattern.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. extractpattern.py:  Call extractptrn module
	1. python3 extractptrn.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --input-pattern: Indicate an input text file and its path to find and filter a specific sequence pattern (accept reverse complement sequences)
		+ --input-length: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (integer only)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)



### asmstatsunlm (mfs)
- To generate a summary of multi-line fasta statistics for unlimited fasta files. An extended version of the allfastats for multiple fasta files. 
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Unlimited multi-line fasta files.
	+ Output: A summary of multi-line fasta with its diverse statistics (e.g. sequence length, GC content, N50, and more).
  	+ Example file: [all_stat_asm.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/all_stat_asm.fa) and [stat_asm_mlt_unlm.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/stat_asm_mlt_unlm.fa) in the "example_data" folder.

Example usage
```
python3 asmstatsunlm.py --input-seqs test_dna1.fasta test_dna2.fasta test_dna3.fasta test_dna4.fasta test_dna5.fasta --out output_dna.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).
+ Both allfastasts and asmstatsunlm modules take the same multi-line FASTA files. Allfastats generates a summary of assembly statistics from a single FASTA file, while asmstatsunlm does so for multiple FASTA files, useful for genome or transcriptome analysis.

- Parameter explanation
	1. python 3: Call python 3
	1. asmstatsunlm.py:  Call asmstatsunlm module
	1. python3 asmstatsunlm.py --help: Check help menu
		+ --input-seqs: Indicate an input multi-line fasta file and its path (Separated by a space for each new fasta file)
		+ --out: Indicate an output text file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### concatenate (ccf)
- To make a concatenated fasta file for unlimited fasta files. 
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multiple fasta files with the same prefix IDs.
	+ Output: A single-line fasta with concatenated IDs and their corresponding sequences.
  	+ Example file: [concat_seq1.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/concat_seq1.fa), [concat_seq2.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/concat_seq2.fa), [concat_seq3.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/concat_seq3.fa), [concat_seq4.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/concat_seq4.fa), and [concat_seq1.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/concat_seq5.fa) in the "example_data" folder.

Example usage
```
python concatenate.py --input-seqs test_dna1.fasta test_dna2.fasta test_dna3.fasta test_dna4.fasta test_dna5.fasta --out output_concat.fasta --t 1 --mem 2
```
+ For optimal use of this module, ensure that all input FASTA files have matching prefix IDs and headers and are formatted as single-line FASTA before concatenation.
+ Before using the concatenate module, ensure your FASTA files are in single-line format by using the renameid and multi2singleline modules, even though the pipeline automatically converts multi-line FASTA files.

- Parameter explanation
	1. python 3: Call python 3
	1. concatenate.py:  Call concatenate module
	1. python3 concatenate.py --help: Check help menu
		+ --input-seqs: Indicate input multiple single-line fasta files and their path (Separated by a space for each new fasta file)
		+ --out: Indicate a concatenated output fasta file and its path
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)
 
### translated (tls)
- To find the translated sequences as a protein including the open reading frames (ORFs)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multline fasta files.
	+ Output: A single-line fasta with translation and selected ORFs (nucleotide and protein sequences).
	+ Example file: [trnsldna.fa](https://github.com/OZTaekOppa/FastQHandler/blob/main/example_data/trnsldna.fa) in the "example_data" folder.

Example usage
```
python translatedna.py --input-seq test_dna.fasta --out output/folder --t 1 --mem 2
```
+ For optimal use of this module, ensure that all input FASTA files have matching prefix IDs and headers and are formatted as single-line FASTA before concatenation.
+ Before using the translatedna module, ensure your FASTA files are in single-line format by using the renameid and multi2singleline modules, even though the pipeline automatically converts multi-line FASTA files.

- Parameter explanation
	1. python 3: Call python 3
	1. translatedna.py:  Call translatedna module
	1. python3 translatedna.py --help: Check help menu
		+ --input-seq: Indicate input single-line fasta files and their path
		+ --out: Indicate a translated output fasta file and its path (a total of four fasta files for both nucleotide and protein sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


## FAQ

We encourage users to use the [Issues](https://github.com/OZTaekOppa/FastQHandler/issues).


## WIKI PAGE

Please see the GitHub page.


## AUTHORS

**Hyungtaek Jung** and the [**National Centre for Indigenous Genomics**](https://ncig.anu.edu.au/).


## COPYRIGHT

The full **FastQHandler** is distributed under the MIT license. 
