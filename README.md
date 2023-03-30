# MCScanx-python-wrapper
This is a python wrapper for running MCScanx

[MCScanX](https://github.com/wyp1125/MCScanX) is an amazing tool for identifying and visualyzing syntenic blocks between two or more genomes.
Unfortunately, the first time I ran it I had a lot of issues, and the manual could be better. So, I decided to write this python wrapper.

### Before running the wrapper, the following tools must be installed and available in the path:

- ```blastp``` (and makeblastdb)
- ```python 3``` or higher (tested on 3.8.12)
- ```Biopython```
- ```MCScanX``` directory **must** be in the working directory and installed. Installation is very simple following the developer's instructions.
- ```The wrapper``` installation is not needed, just download the code and is ready to be used.

All of the dependencies can be easily installed using developer's instructions or ```conda```.

### Input data:
- Genome references must be in fasta format. Chromosome names **must** contain 2 letters and one number (i.e ch1, ch2, ch3, etc)
- Genome annotations must be in gff/gff3 format (gene id must be indicated in the 8th column and with the following format: ID=). The chromosome name must match the gff.

Now, you're ready to run it!

### Ussage of the wrapper:
```
python macscanxWrapper.py -h 
usage: macscanxWrapper.py [-h] -f1 specie1.fasta -g1 specie1.gff -f2 file2.fasta -g2 specie2.gff -t threads(int) -p prefix

	Python wrapper for MCScanX
---------------------
Tested on python 3.8.12

optional arguments:
  -h, --help            show this help message and exit

Mandatory arguments:
  -f1, --fasta1         Fasta file for specie1 (Genomic sequence)
  -g1, --gff11          GFF/GFF3 file for specie1
  -f2, --fasta2         Fasta file for specie2 (Genomic sequence)
  -g2, --gff2           GFF/GFF3  file for specie2
  -t, --blastp_threads  Number of threads (CPUs) to use in the BLAST search
  -p, --prefix          Job prefix, it will be used for MCScanX results
```

## The pipeline:
![pythonWrapper](https://user-images.githubusercontent.com/45425927/228645462-4ad97d34-0448-478a-815e-522ef74b485a.png)

### Step 1. Checking dependencies

The first step is to check all dependencies are in the correct path and/or running.

### Step 2. Loading/parsing input data

- By default, only **gene** features will be extracted.
- By default the standard code (table 1) is used for protein translation.
- Pseudogenes (missing any start codon) will be ignored.
- MCScanX needs a bed-like file containing the annotation (with the following format: **chr**\t**geneID**\t**Start**\t**End**), it is created and named as **prefix.gff** (please notice it is not a conventional gff file, but is used for visualization i.e using [SynVisio](https://synvisio.github.io/#/)).

### Step 3. Running blastp and MCScanX

- Proteins from species1 and species2 will be merged into a single file. If you want to compare >2 species, there is no need to duplicate them in the input data, they can be present either in species1 or species2.
- By default, blastp will only keep the 5 matches for every sequence in the query.
- The number of threads for running blastp must be specified (reccomended >= 15).
- MCScanX will run in its own directory.

### Step 4. Output files

- Output files from MCScanX will be placed in the working directory.

### Running the wrapper using test files:

```
python macscanxWrapper.py -f1 g3_test.fasta -g1 g3_test.gff3 -f2 3688_test.fasta -g2 3688_test.gff3 -t 15 -p test_final

Checking dependencies... 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
./MCScanX-master directory is ready...
./MCScanX is ready...
blastp command is ready...
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Loading and parsing input files...

Valid sequences (proteins) in g3_test.fasta: 11195
Valid Sequences (proteins) in 3688_test.fasta: 10999
Annotation file ready: test_final.gff

Generate blastdb:
makeblastdb -in species_merged.faa -dbtype prot
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 


Building a new DB, current time: 03/29/2023 12:52:37
New DB name:   /Users/Documents/python_wrapperMCScanX/species_merged.faa
New DB title:  species_merged.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 22194 sequences in 0.444309 seconds.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Running Blastp:
blastp -query species_merged.faa -db species_merged.faa -outfmt 6 -out test_final.blast -max_hsps 5 -max_target_seqs 5 -num_threads 15
Files for MCScanX are ready: test_final.blast, test_final.gff

MCScanX is running ... 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
Reading BLAST file and pre-processing
Generating BLAST list
66695 matches imported (40987 discarded)
3 pairwise comparisons
249 alignments generated
Pairwise collinear blocks written to test_final.collinearity [2.978 seconds elapsed]
Tandem pairs written to test_final.tandem
Writing multiple syntenic blocks to HTML files
tv1.html
vt1.html
Done! [0.582 seconds elapsed]
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
MCScanX finished succesfully ... 
```
### Output files 

- test_final.blast
- test_final.collinearity
- test_final.gff
- test_final.html **(directory)**
- test_final.tandem



