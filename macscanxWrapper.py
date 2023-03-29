#!/usr/bin/env python
#-*- coding: UTF-8 -*-
from __future__ import division
import sys
import argparse
import operator
import re
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
import shutil

program = 'macscanxWrapper.py'
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description="\n\n\tPython wrapper for MCScanX\n---------------------\nTested on python 3.8.12\n\n")
requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-f1", "--fasta1", dest='fa_1', required=True, type=str, help='Fasta file for specie1 (Genomic sequence)')
requiredNamed.add_argument("-g1", "--gff11", dest='gff_1', required=True, type=str, help='GFF/GFF3 file for specie1')
requiredNamed.add_argument("-f2", "--fasta2", dest='fa_2', required=True, type=str, help='Fasta file for specie2 (Genomic sequence)')
requiredNamed.add_argument("-g2", "--gff2", dest='gff_2', required=True, type=str, help='GFF/GFF3 file for specie2')
requiredNamed.add_argument("-t", "--blastp_threads", dest='bpthreads', required=True, type=int, help='Number of threads (CPUs) to use in the BLAST search')
requiredNamed.add_argument("-p", "--prefix", dest='wprefix', required=True, type=str, help='Job prefix, it will be used for MCScanX results')
#requiredNamed.add_argument("-b", "--bed", dest='ann_bed', required=True, type=str, help='Bed-like file containing annotations (chrom\tgeneID\tstart\tend), and ending by .gff')
args = parser.parse_args()

##########################################################################################
##########################################################################################
#Bed-like name must match the prefix;

print('\nChecking dependencies... ')
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
#MCScanX directory must exist in the same path

if os.path.exists('./MCScanX-master'):
	print('./MCScanX-master directory is ready...')
else:
	print('FATAL ERROR: MCScanX directory is not in the working directory!  ... ')
	sys.exit()

if shutil.which('./MCScanX-master/MCScanX') is not None:
	print('./MCScanX is ready...')
else:
	print('FATAL ERROR: MCScanX is not correctly installed!  ... ')
	sys.exit()

if shutil.which('blastp') is not None:
	print('blastp command is ready...')
else:
	print('FATAL ERROR: blastp is not correctly installed!  ... ')
	sys.exit()
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')

##########################################################################################
##########################################################################################


##########################################################################################
#Step 1; Extract protein sequences from annotation

species_list = [args.fa_1, args.gff_1, args.fa_2, args.gff_2]

#Reverse complement sequence, function
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def reverse_complement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

#nucleotide to protein function
def nucltoprot(gff_file, fastafile):
	seq = ''
	fasta_seqs = {}
	with open(fastafile, 'r') as fa_in:
		for line in fa_in:
			line = line.rstrip('\n')
			if line.startswith('>'):
				if seq:
					fasta_seqs[name] = seq
					seq = ''
				name = line[1:].strip()
			else:
				seq = seq + line
		fasta_seqs[name] = seq
	faa_dict = {}
	with open(gff_file, 'r') as gff_in:
		for line in gff_in:
			line = line.rstrip('\n')
			col = line.split('\t')
			if col[2] == 'gene':
				id_gene = col[8].split('=')[1]
				if ';' in id_gene:
					id_gene = col[8].split('=')[1].split(';')[0]
				else:			
					id_gene = id_gene
				seq_fasta = fasta_seqs[col[0]]
				if col[6] == '+':
					#faa_dict[id_gene] = Seq(seq_fasta[int(col[3])-1 :int(col[4])]).translate(to_stop=True, gap='N')
					try:
						faa_dict[id_gene] = Seq(max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',(seq_fasta[int(col[3])-1 :int(col[4])])), key = len)).translate(cds=True, gap='N')
					except ValueError:
						continue
				elif col[6] == '-':
					try:
						#faa_dict[id_gene] = Seq(reverse_complement(seq_fasta[int(col[3])-1 :int(col[4])])).translate(to_stop=True, gap='N')
						faa_dict[id_gene] = Seq(max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',reverse_complement(seq_fasta[int(col[3])-1 :int(col[4])])), key = len)).translate(cds=True, gap='N')
					except ValueError:
						continue
				else:
					try:
						#faa_dict[id_gene] = Seq(seq_fasta[int(col[3])-1 :int(col[4])]).translate(to_stop=True, gap='N')
						faa_dict[id_gene] = Seq(max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',seq_fasta[int(col[3])-1 :int(col[4])]), key = len)).translate(cds=True, gap='N')
					except ValueError:
						continue
	return faa_dict

ct_sp1 = 0
ct_sp2 = 0
print('\nLoading and parsing input files...')
with open(args.fa_1.split('.')[0] + '_prot.faa', 'w') as ouprot1:
	for k, v in nucltoprot(species_list[1], species_list[0]).items():
		ct_sp1 += 1
		ouprot1.write('>' + str(k) + '\n' + str(v) + '\n')
with open(args.fa_2.split('.')[0] + '_prot.faa', 'w') as ouprot_2:
	for k, v in nucltoprot(species_list[3], species_list[2]).items():
		ct_sp2 += 1
		ouprot_2.write('>' + str(k) + '\n' + str(v) + '\n')

print('\nValid sequences (proteins) in ' + str(args.fa_1) + ': ' + str(ct_sp1))
print('Valid Sequences (proteins) in '  + str(args.fa_2) + ': ' + str(ct_sp2))

filenames = [args.fa_1.split('.')[0] + '_prot.faa', args.fa_2.split('.')[0] + '_prot.faa']
with open('species_merged.faa', 'w') as outMspecies:
	for fname in filenames:
		with open(fname) as infile:
			for line in infile:
				outMspecies.write(line)

##########################################################################################
#Create annotation files from gff, it is a bed-like format chr\tgene_id\tstart\tend
#discarded coordinates (if any) are not discarded for the annotation

def gff_to_mcscanformat(gff_file):
	gff_dict = {}
	with open(gff_file, 'r') as gff_in:
		for line in gff_in:
			line = line.rstrip('\n')
			col = line.split('\t')
			if col[2] == 'gene':
				id_gene = col[8].split('=')[1]
				if ';' in id_gene:
					id_gene = col[8].split('=')[1].split(';')[0]
				else:			
					id_gene = id_gene
				gff_dict[id_gene] = [str(col[0]), str(col[3]), str(col[4])]
	return gff_dict

species_list = [args.fa_1, args.gff_1, args.fa_2, args.gff_2]

with open(args.wprefix + '.gff', 'w') as outanot:
	for k, v in gff_to_mcscanformat(species_list[1]).items():
		outanot.write('\t'.join([str(v[0]), str(k), str(v[1]), str(v[2])])  + '\n')
	for k, v in gff_to_mcscanformat(species_list[3]).items():
		outanot.write('\t'.join([str(v[0]), str(k), str(v[1]), str(v[2])])  + '\n')

annot_file = args.wprefix + '.gff'

print('Annotation file ready: ' + str(annot_file))

##########################################################################################
#Step 2; Create the blastp database

command = "makeblastdb -in species_merged.faa -dbtype prot"
print('\nGenerate blastdb:')
print("makeblastdb -in species_merged.faa -dbtype prot")
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
subprocess.call(command, shell=True)

##########################################################################################
#Step 3; Running blastp
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')

print('\nRunning Blastp:')
command_2 = 'blastp -query species_merged.faa -db species_merged.faa -outfmt 6 -out '+ str(args.wprefix) + '.blast -max_hsps 5 -max_target_seqs 5 -num_threads ' + str(args.bpthreads)
print(command_2)
subprocess.call(command_2, shell=True)

##########################################################################################
#Running MCScanX

print('Files for MCScanX are ready: ' + str(args.wprefix) + '.blast, ' + annot_file)
subprocess.call("mv " + str(args.wprefix) + '.blast ./MCScanX-master', shell=True)
subprocess.call("mv " + annot_file + ' ./MCScanX-master', shell=True)
os.chdir("./MCScanX-master")
print('\nMCScanX is running ... ')
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
subprocess.call("./MCScanX " + str(args.wprefix), shell=True)
subprocess.call('mv ./' + str(args.wprefix) + '.* ./../', shell=True)
print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')

print('MCScanX finished succesfully ... \n\n')


