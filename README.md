# MCScanx-python-wrapper
This is a python wrapper for running MCScanx

[MCScanX](https://github.com/wyp1125/MCScanX) is an amazing tool for identifying and visualyzing syntenic blocks between two or more genomes.
Unfortunately, the first time I ran it I had a lot of issues, and the manual is not the best. So I decided to create this python wrapper.


### Before running the wrapper, the following tools must be installed and available in the path:

```
blastp (and makeblastdb)
python 3 or higher (tested on 3.8.12)
Biopython
MCScanX directory MUST be in the working directory and installed. Installation is very simple following the developer's instructions.
The wrapper: installation is not needed, just download the code and is ready to be used.
```
### Input data:
1.- Genome references must be in fasta format. Chromosome names ***must*** contain 2 letters and one number (i.e ch1)
2.- Genome annotations must be in gff/gff3 format (gene id must be indicated in the 8th column and with the following format: ID=). The chromosome name must match the gff.

Now you're ready to run it!

## The pipeline:
![pythonWrapper](https://user-images.githubusercontent.com/45425927/228645462-4ad97d34-0448-478a-815e-522ef74b485a.png)

### Step 1. Checking dependencies

The first step is to check all dependencies are in the correct path and/or running.

### Step 2. Loading/parsing input data

- By default, only ***gene*** features will be extracted.
- For protein translation the standard code (table 1) will be used.
- Pseudogenes (missing any start codon) will be ignored.
- MCScanX needs a bed-like file containing the annotation (chr\tgeneID\tStart\tEnd), it is created and named as ***prefix.gff***.

### Step 3. Running blastp and MCScanX

- By default, blastp will only keep the 5 matches for every sequence in the query.
- The number of threads for running blastp must be specified (reccomended >= 15).
- MCScanX will run in its own directory.

### Step 4. Output files

- Output files from MCScanX will be placed in the working directory.







