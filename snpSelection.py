#!/usr/bin/env python
##source from Hui Yang.

import os
import sys
import pysam
import twobitreader
import math
import subprocess
import pandas as pd
from scipy.stats import binom_test
from scipy.stats import fisher_exact

def split_bam_file(bamF,snpF,strain1,strain2,mix):
	snp = {}
	for line in open(snpF,'r'):
		if line.startswith("#"):
			continue
		line = line.strip().split('\t')
		if not line[0].startswith("chr"):
			line[0]=f"chr{line[0]}"
		if line[0] not in snp:
				snp[line[0]] = {}
		snp[line[0]][int(line[1])] = [line[3],line[4]]

	votes = {}
	match, unmatch = 0,0
	infile = pysam.Samfile(bamF,'rb')
	for line in infile:
		if line.qname not in votes:
				votes[line.qname] = [0,0]
		# print(line.rname)
		chr = line.reference_name
		start = line.pos + 1
		if line.cigar:
			read = ''
			score = ''
			index = 0
			for cg in line.cigar:
				if cg[0] == 0:
					read += line.seq[index:(index+cg[1])]
					score += line.qual[index:(index+cg[1])]
					index += cg[1]
				elif cg[0] == 2:
					read += 'N'*cg[1]
					score += '!'*cg[1]
					index += cg[1]
				else:
					index += cg[1]
			for pos in range(0,len(read)):
				if chr in snp and start+pos in snp[chr] and ord(score[pos])-33 >= 30:
					if read[pos] in snp[chr][start+pos][0].split(';'): # A;T
						votes[line.qname][0] += 1
						match += 1
					elif read[pos] in snp[chr][start+pos][1].split(';'): # C;G
						votes[line.qname][1] += 1
						match += 1
					else:
						unmatch += 1

	outf1 = pysam.Samfile(strain1,'wb',template=infile)
	outf2 = pysam.Samfile(strain2,'wb',template=infile)
	outf3 = pysam.Samfile(mix,'wb',template=infile)
	maternal, paternal, mixed = 0, 0, 0
	for line in pysam.Samfile(bamF,'rb'):
		if line.qname in votes:
			if sum(votes[line.qname]) > 0:
				if votes[line.qname][0]*1.0/sum(votes[line.qname]) >= 2.0/3:
					maternal += 1
					outf1.write(line)
				elif votes[line.qname][0]*1.0/sum(votes[line.qname]) <= 1.0/3:
					paternal += 1
					outf2.write(line)
				else:
					mixed += 1
					outf3.write(line)
	outf1.close()
	outf2.close()
	outf3.close()
	infile.close()

	print(match,unmatch)
	print(maternal, paternal, mixed)


def split_bam_file_methyl(bamF,snpF,strain1,strain2):
	snp = {}
	for line in open(snpF,'r'):
		line = line.strip().split('\t')
		if line[0] not in snp:
				snp[line[0]] = {}
		snp[line[0]][int(line[1])] = [line[3],line[4]]

	votes = {}
	match, unmatch = 0,0
	infile = pysam.Samfile(bamF,'rb')
	for line in infile:
		if line.qname not in votes:
				votes[line.qname] = [0,0]
		# print(line.rname)
		chr = line.reference_name
		start = line.pos + 1
		if line.cigar:
			read = ''
			score = ''
			index = 0
			for cg in line.cigar:
				if cg[0] == 0:
					read += line.seq[index:(index+cg[1])]
					score += line.qual[index:(index+cg[1])]
					index += cg[1]
				elif cg[0] == 2:
					read += 'N'*cg[1]
					score += '!'*cg[1]
					index += cg[1]
				else:
					index += cg[1]
			for pos in range(0,len(read)):
				if chr in snp and start+pos in snp[chr] and ord(score[pos])-33 >= 30:
					if read[pos] in snp[chr][start+pos][0].split(';'):
						votes[line.qname][0] += 1
						match += 1
					elif read[pos] in snp[chr][start+pos][1].split(';'):
						votes[line.qname][1] += 1
						match += 1
					else:
						unmatch += 1

	outf1 = pysam.Samfile(bamF[:-4]+'_'+strain1+'.bam','wb',template=infile)
	outf2 = pysam.Samfile(bamF[:-4]+'_'+strain2+'.bam','wb',template=infile)
	outf3 = pysam.Samfile(bamF[:-4]+'_mixed.bam','wb',template=infile)
	maternal, paternal, mixed = 0, 0, 0
	for line in pysam.Samfile(bamF,'rb'):
		if line.qname in votes:
			if sum(votes[line.qname]) > 0:
				if votes[line.qname][0]*1.0/sum(votes[line.qname]) >= 2.0/3:
					maternal += 1
					outf1.write(line)
				elif votes[line.qname][0]*1.0/sum(votes[line.qname]) <= 1.0/3:
					paternal += 1
					outf2.write(line)
				else:
					mixed += 1
					outf3.write(line)
	outf1.close()
	outf2.close()
	outf3.close()
	infile.close()

	print(match,unmatch)
	print(maternal, paternal, mixed)


def allelic_MPBN(infile, outfile, c1, c2, valid_count = 10):
	df = pd.read_table(infile, sep ='\t', header = None, skiprows = 0)
	# df.columns = ['#Chr','Start','End','RefSeq','GeneSymbol','Strand','txStart','txEnd','Score','exonCount','exonStarts','exonEnds','BP.E65_EPI.Mat','BP.E65_EPI.Pat','PB.E65_EPI.Mat','PB.E65_EPI.Pat','BP.E65_VE.Mat','BP.E65_VE.Pat','PB.E65_VE.Mat','PB.E65_VE.Pat','BP.E65_EXE.Mat','BP.E65_EXE.Pat','PB.E65_EXE.Mat','PB.E65_EXE.Pat','BP.E95_Placenta.Mat','BP.E95_Placenta.Pat','PB.E95_Placenta.Mat','PB.E95_Placenta.Pat']
	outf = open(outfile, 'w')
	print('RefSeq\tGeneSymbol\tMat\tPat\tp-value\tAS\tclass', file=outf)

	with open(infile, 'r') as f:
		# next(f)
		# next(f)
		for line in f:
			line = line.strip().split('\t')
			if int(line[c1]) + int(line[c2]) < valid_count:
				# continue
				print('\t'.join([line[3],line[4],line[c1],line[c2],'-','-','N']), file=outf)
			else:
				# pBase = float(df.iloc[:,c1].sum())/(float(df.iloc[:,c1].sum()) + float(df.iloc[:,c2].sum()))
				# pvalue = binom_test(int(line[c1]),int(line[c1]) + int(line[c2]), pBase)+1E-300
				pvalue = binom_test(int(line[c1]),int(line[c1]) + int(line[c2]), 0.5)+1E-300
				if pvalue <= 0.001:
					if int(line[c1]) > int(line[c2]):
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'M']), file=outf)
					else:
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'P']), file=outf)
				else:
					if int(line[c1]) > int(line[c2]):
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'B']), file=outf)
					else:
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'B']), file=outf)
		outf.close()

def mkdirs(path):
	path=os.path.abspath(path)
	if os.path.isfile(path):
		p=os.path.split(path)[0]
	if not os.path.exists(path):
		os.makedirs(path)
	return True


def main():
	selection = sys.argv[1]
	for i in sys.argv[3:7]:
		paths=i.split("/")
		if len(paths) ==1:
			continue
		parent_path=os.path.join(*paths[:-1])
		print(parent_path)
		mkdirs(parent_path)
	if selection == "split":
		split_bam_file(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]) # bamF, snpF, strain1, strain2, mixname
	if selection == "splitMethyl":
		split_bam_file_methyl(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) # bamF, snpF, strain1, strain2
	if selection == "allele_negLogPValue":
		allelic_MPBN(sys.argv[2],sys.argv[3], int(sys.argv[4]), int(sys.argv[5])) #infile, outfile, c1, c2, valid_count


main()
