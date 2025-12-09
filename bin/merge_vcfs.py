#!/usr/bin/python3

# MIT License

# Copyright (c) 2020 Jana Ebler

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import sys
import argparse
import re
import pyfaidx
from collections import defaultdict

class Variant:

	def __init__(self, samples, start, ref_allele, alt_alleles, genotypes, ploidy, var_ids):
		assert len(samples) == len(genotypes)
		assert len(alt_alleles) == len(var_ids)
		self._samples = samples
		self._start = start
		self._end = start + len(ref_allele)
		self._ref_allele = ref_allele
		self._alt_alleles = alt_alleles
		self._genotypes = genotypes
		self._ploidy = ploidy
		self._id = var_ids

	def get_samples(self):
		return self._samples

	def get_start(self):
		return self._start

	def get_end(self):
		return self._end

	def get_ref(self):
		return self._ref_allele

	def get_alt(self):
		return self._alt_alleles

	def get_ploidy(self):
		return self._ploidy

	def get_id(self, allele_index):
		if allele_index == 0:
			assert(False)
		else:
			return self._id[allele_index-1]

	def get_ids(self):
		return self._id

	def get_allele(self, index):
		if index == 0:
			return self._ref_allele
		else:
			return self._alt_alleles[index - 1]

	def get_genotype_of(self, sample):
		if sample in self._samples:
			return self._genotypes[self._samples.index(sample)]
		else:
			return Genotype([0] * self._ploidy, True)

	def is_absent(self):
		alleles = []
		for g in self._genotypes:
			for a in g.get_alleles():
				alleles.append(a)
		return all(x == 0 for x in alleles)

	def __eq__(self, other):
		if self._start != other._start:
			return False
		if self._end != other._end:
			return False
		if self._samples != other._samples:
			return False
		if self._ref_allele != other._ref_allele:
			return False
		if self._alt_alleles != other._alt_alleles:
			return False
		if self._genotypes != other._genotypes:
			return False
		if self._ploidy != other._ploidy:
			return False
		if self._id != other._id:
			return False
		return True

	def __lt__(self, other):
		return self._start < other._start

	def __le__(self, other):
		return self._start <= other._start

	def __gt__(self,other):
		return self._start > other._start

	def __ge__(self, other):
		return self._start >= other._start

	def __str__(self):
		sample_str = '[' + ','.join([str(s) for s in self._samples]) + ']'
		start_str = str(self._start)
		end_str = str(self._end)
		alt_str = ','.join(self._alt_alleles)
		genotype_str = '[' + ','.join([str(g) for g in self._genotypes]) + ']'
		id_str = '[' + ','.join(self._id) + ']'
		return 'Variant(' + ','.join([sample_str, start_str, end_str, self._ref_allele, alt_str, genotype_str, id_str]) + ')'


class Allele:
	"""
	Represents an alternative allele.
	"""

	def __init__(self, allele, chrom, start, end, id):
		"""
		allele: str
			alternative allele sequence
		chrom: str
			chromosome this allele is located on
		start, end: int
			starting and end position relative
			to reference sequence
		id: str
			allele ID
		"""	
		self._allele = allele
		self._start = start
		self._id = id
		self._end = end
		self._chrom = chrom

	def chrom(self):
		return self._chrom

	def start(self):
		return self._start

	def end(self):
		return self._end

	def allele(self):
		return self._allele

	def id(self):
		return self._id

	def __str__(self):
		return 'Allele(' + ','.join([str(self._allele), self._chrom, str(self._start), str(self._end), str(self._id)]) + ')'


def insert_allele_sequence(sequence, allele, offset):
	"""
	Given a sequence (reference) insert the given allele.
	"""
	start = allele.start() - offset
	end = allele.end() - offset
	included = sequence[0:start] + allele.allele() + sequence[end:]
	offset -= len(included) - len(sequence)
	return included, offset


class Haplotype:
	"""
	Represents a haplotype.
	"""

	def __init__(self, alleles):
		"""
		alleles: list(int)
			list of alleles along the haplotype
		"""
		self._alleles = tuple(alleles)

	def is_undefined(self):
		"""
		Checks whether the haplotype is undefined.
		This is the case if at least one of its alleles
		is missing (-1).
		"""
		return -1 in self._alleles

	def is_reference(self):
		"""
		Checks whether the haplotype contains only
		reference alleles.
		"""
		return all([h == 0 for h in self._alleles])

	def __eq__(self, other):
		return self._alleles == other._alleles

	def __lt__(self, other):
		return self._alleles < other._alleles

	def __hash__(self):
		return hash(self._alleles)

	def get_alleles(self):
		return self._alleles
		
	

class Genotype:
	"""
	Represents a genotype.
	"""

	def __init__(self, alleles, is_phased):
		"""
		alleles: list
			alleles of the genotype
		is_phased: bool
			whether the genotype is phased.
		"""
		self._alleles = alleles
		self._phased = is_phased

	def __str__(self):
		allele_str = []
		for allele in self._alleles:
			if allele == -1:
				allele_str.append('.')
			else:
				allele_str.append(allele)
		if self._phased:
			return '|'.join([str(i) for i in allele_str])
		else:
			return '/'.join([str(i) for i in allele_str])

	def get_alleles(self):
		return self._alleles

	def get_ploidy(self):
		return len(self._alleles)

	def is_phased(self):
		return self._phased

	def __eq__(self, other):
		if self._phased != other._phased:
			return False
		if self._phased:
			return self._alleles == other._alleles
		else:
			return sorted(self._alleles) == sorted(other._alleles)

	def is_hom_ref(self):
		return all([int(a) == 0 for a in self._alleles])



def genotype_from_string(gt_string):
	"""
	Given genotype string from VCF, convert into a Genotype
	object.
	"""
	is_phased = False
	alleles = []
	if '|' in gt_string:
		is_phased = True
		alleles = []
		for allele in gt_string.split('|'):
				if allele != '.':
					alleles.append(int(allele))
				else:
					alleles.append(-1)
	elif '/' in gt_string:
		is_phased = False
		for allele in gt_string.split('/'):
				if allele != '.':
					alleles.append(int(allele))
				else:
					alleles.append(-1)
	else:
		is_phased = True
		alleles = [int(gt_string)] if gt_string != '.' else [-1]
	return Genotype(alleles, is_phased)


def genotypes_from_list(alleles, ploidy):
	"""
	Given a list of alleles and a ploidy,
	construct the genotypes.
	"""
	genotypes = []
	n_alleles = len(alleles)
	assert n_alleles % ploidy == 0
	for i in range(0, n_alleles, ploidy):
		genotypes.append('|'.join([str(a) for a in alleles[i:i + ploidy]]))
	return genotypes


class HaplotypeTable:
	"""
	Represents columns of Haplotypes.
	
	Attributes:
		haplotypes: maps a column index to the alleles along
			haplotypes of this column
		ploidy: ploidy of the data
		chrom: chromosome
		alleles: maps row index to Alleles
	"""

	def __init__(self, ploidy):
		"""
		ploidy: int
			the ploidy of the samples
		"""
		self._ploidy = ploidy
		self._chrom = None
		self._haplotypes = defaultdict(list)
		self._alleles = defaultdict(list)
		self._start = float("inf")
		self._end = 0

	def is_empty(self):
		"""
		Check whether table is empty.
		"""
		return not self._alleles

	def insert_allele_index(self, col_index, allele):
		"""
		Add an allele to the haplotype in column col_index.
		"""
		self._haplotypes[col_index].append(allele)

	def insert_allele_sequence(self, row_index, allele):
		"""
		Insert an allele for row index row_index.
		"""
		self._alleles[row_index].append(allele)
		self._start = min(self._start, allele.start())
		self._end = max(self._end, allele.end())
		if self._chrom is not None:
			assert self._chrom == allele.chrom()
		self._chrom = allele.chrom()

	def parse(self, row_index, line, id_in_info=False):
		"""
		Parse line from VCF file and store sample columns.
		"""
		total_ids = set([])
		fields = line.split()
		# parse allele IDs either from ID column or ID field in info column (ids need to
		# be separated by a comma
		ids = fields[2].split(',')
		
		ids = fields[2].split(',')
		if id_in_info:
			info = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
			assert 'ID' in info
			ids = info['ID'].split(',')
		alt_alleles = fields[4].split(',')
		if len(alt_alleles) != len(ids):
			raise Exception('An ID needs to be provided for each individual ID, either in ID colum or info column.')

		# store the alleles (and their IDs) in HaplotypeTable
		for id,allele in zip(ids,alt_alleles):
			start = int(fields[1])
			end = start + len(fields[3])
			allele_obj = Allele(allele, fields[0], start, end, id)
			self.insert_allele_sequence(row_index, allele_obj)
			total_ids.add(id)
			

		# store the genotype alleles of each sample in HaplotypeTable
		column_index = 0
		for genotype_string in fields[9:]:
			format_field = fields[8].split(':')
			if not 'GT' in format_field:
				raise Exception('Genotype information missing from VCF file.')
			gt_index = format_field.index('GT')
			genotype = genotype_from_string(genotype_string.split(':')[gt_index])
			if not genotype.is_phased():
				raise Exception('Line: ' + line + ' contains an unphased genotype.')
			for allele in genotype.get_alleles():
				self.insert_allele_index(column_index, allele)
				column_index += 1
		return len(total_ids)


	def merge(self, reference):
		"""
		Merge all rows into a single variant combining all alleles and
		print the result as a VCF record.
		"""

		written_ids = set([])
		
		# no variants present, return
		if self.is_empty():
			return None, 0
		
		# reference sequence of this variant cluster
		ref_allele = reference[self._chrom][self._start-1 : self._end-1]
		# invert map
		hap_to_column = defaultdict(list)
		for k, v in self._haplotypes.items():
			hap_to_column[Haplotype(v)].append(k)

		# construct the haplotype sequences covered by samples in this region
		hap_to_sequence = {}
		hap_to_id = {}
		for hap in hap_to_column.keys():
			prev_end = 0
			undefined_allele = hap.is_undefined()
			if not undefined_allele:
				sequence = ref_allele
				id_sequence = []
				offset = self._start
				for row_index, allele in enumerate(hap.get_alleles()):
					assert allele != -1
					if allele > 0:
						start = self._alleles[row_index][allele-1].start()
						if prev_end > start:
							sys.stderr.write('Two overlapping variants at same haplotype at ' + self._chrom + ':' + str(self._start) + ', set allele to missing.\n')
							undefined_allele = True
							break
						sequence, offset = insert_allele_sequence(sequence, self._alleles[row_index][allele-1], offset)
						id_sequence.append(self._alleles[row_index][allele-1].id())
						prev_end = self._alleles[row_index][allele-1].end()
				if not undefined_allele:
					hap_to_sequence[hap] = sequence
					hap_to_id[hap] = ':'.join(id_sequence)
		
		# construct VCF line representing the merged variant
		genotypes = ['.'] * len(self._haplotypes)
		alt_alleles = []
		ids = []
		allele_index = 1
		for haplotype in hap_to_sequence.keys():
			allele = allele_index
			if haplotype.is_reference():
				allele = 0
			elif hap_to_sequence[haplotype] in alt_alleles:
				allele = alt_alleles.index(hap_to_sequence[haplotype]) + 1
				if not hap_to_id[haplotype] in ids:
					sys.stderr.write('Different allele combinations lead to same sequence at ' + self._chrom + ':' + str(self._start) + '.\n')
			else:
				allele = allele_index
				alt_alleles.append(hap_to_sequence[haplotype])
				ids.append(hap_to_id[haplotype])
				for id in hap_to_id[haplotype].split(':'):
					written_ids.add(id)			
				allele_index += 1

			for col in hap_to_column[haplotype]:
				genotypes[col] = allele
		
		# check if no alternative sequences left (e.g. because all haplotypes contained unknown alleles)
		if len(alt_alleles) == 0:
			assert not written_ids
			return None, 0
		vcf_alt = ','.join(alt_alleles)
		
		# check if there are any alleles with N characters
		if any(c not in 'CAGTcagt,' for c in ref_allele) or any(c not in 'CAGTcagt,' for c in vcf_alt):
			return None, 0
		vcf_genotypes = genotypes_from_list(genotypes, self._ploidy)
		vcf_line = [	self._chrom, # CHROM
				str(self._start), # POS
				'.', # ID
				ref_allele, # REF
				','.join(alt_alleles), # ALT
				'.', # QUAL
				'PASS', # FILTER
				'ID=' + ','.join(ids), # INFO 
				'GT', # FORMAT
				'\t'.join([str(g) for g in vcf_genotypes]) # sample columns
			]
		return '\t'.join(vcf_line), len(written_ids)


def parse_line(samples, line, id_in_info=False):
	variants = []
	fields = line.split()
	chrom = fields[0]
	start = int(fields[1])
	ref_allele = fields[3]
	alt_alleles = fields[4].split(',')
	ids = fields[2].split(';')
	if id_in_info:
		# pav vcfs have ids in INFO column...
		info = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
		assert 'ID' in info
		ids = info['ID'].split(',')
	for i,s in enumerate(samples):
		# get samples genotype
		genotype = genotype_from_string(fields[9 + i])
		if not genotype.is_phased():
			raise Exception('Line: ' + line + ' contains an unphased genotype.')
		variants.append(Variant([s], start, ref_allele, alt_alleles, [genotype], genotype.get_ploidy(), ids))
	return chrom, variants
		

def combine_haplotypes(variants):
	alleles = []
	for v in variants:
		samples = v.get_samples()
		assert len(samples) == 1
		genotype = v.get_genotype_of(samples[0])
		alleles += genotype.get_alleles()
	return Genotype(alleles, True)


def print_header(samples):
	print('##fileformat=VCFv4.2')
	print('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">')
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
	print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))


def run_merge(reference, vcf, ploidy, chromosomes=None):
	"""
	Runs variant merging.
	
	Parameters:
	
	reference: str
		name of the reference FASTA file.
	vcf: str
		name of the VCF file.
	ploidy: int
		ploidy of the samples.
	chromosomes:
		only output variants on these chromosomes.
	"""

	# statistics per variant type
	total_input = 0
	total_written = 0

	ref = pyfaidx.Fasta(reference, as_raw=True, sequence_always_upper=True)
	samples = None
	prev_start = float("inf")
	prev_end = 0
	prev_chrom = None
	haplotypes = HaplotypeTable(ploidy)
	row_index = 0
	# read vcf-file
	for line in open(vcf, 'r'):
		if line.startswith('##'):
			continue
		fields = line.split()
		if line.startswith('#'):
			if samples is not None:
				 raise Exception('File ' + filename + ' is not a valid VCF file.')
			if len(fields) < 10:
				raise Exception('File ' + filename + ' does not contain any samples.')
			samples = fields[9:]
			print_header(samples)
			continue
		chromosome = fields[0]
		if chromosomes and (not chromosome in chromosomes):
			continue
		start = int(fields[1])
		ref_allele = fields[3]
		alt_alleles = fields[4].split(',')
		end = start + len(ref_allele)
		if ((start >= prev_end) or (prev_chrom != chromosome)) and not haplotypes.is_empty():
			# print current cluster
			vcf_line, written = haplotypes.merge(ref)
			if vcf_line is not None:
				print(vcf_line)
			total_written += written
			haplotypes = HaplotypeTable(ploidy)
			row_index = 0
		total_input += haplotypes.parse(row_index, line)
		row_index += 1
		prev_end = max(end, prev_end) if (prev_chrom == chromosome) else end
		prev_chrom = chromosome
	if not haplotypes.is_empty():
		vcf_line, written = haplotypes.merge(ref)
		if vcf_line is not None:
			print(vcf_line)
		total_written += written
	
	# print statistics
	sys.stderr.write('Total number of input alleles: ' + str(total_input) + '\n')
	sys.stderr.write('Total number of written alleles: ' + str(total_written) + '\n')


def run_combine_columns(vcf, samples):
	"""
	Combines columns beloning to the same sample.
	
	Parameters:
	
	vcf: str
		name of the input VCF file
	samples: str
		file specifying which columns of the VCF
		belong to the same sample and shall be
		combined
	
	"""

	sample_to_haplotype = defaultdict(list)
	for line in open(samples, 'r'):
		fields = line.split()
		sample_to_haplotype[fields[0]] = fields[1:]
	column_to_index = {}
	with open(vcf, 'r') as outfile:
		nr_samples = None
		for line in outfile:
			if line.startswith('##'):
				print(line[:-1])
				continue
			if line.startswith('#'):
				columns = line[:-1].split('\t')[9:]
				for index,c in enumerate(columns):
					column_to_index[c] = index
				samples = [s for s in sample_to_haplotype.keys()]
				print('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples))	
				continue
			vcf_samples = []
			for s,v in sample_to_haplotype.items():
				vcf_samples.extend(v)
			chromosome, variants = parse_line(vcf_samples,line,False)
			fields = line.split('\t')[:9]
			fields[8] = 'GT'
			for sample,haps in sample_to_haplotype.items():
				fields.append(str(combine_haplotypes([variants[column_to_index[h]] for h in haps])))
			print('\t'.join(fields))




if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='merge_vcfs.py', description=__doc__)
	subparsers = parser.add_subparsers(dest="subparser_name")
	parser_merge = subparsers.add_parser('merge', help='merge variant calls into pangenome graph represented in terms of a multi-sample VCF file.')
	parser_merge.add_argument('-r', metavar='reference', required=True, help='reference sequence (FASTA-format).')
	parser_merge.add_argument('-vcf', metavar='VCF', required=True, help='multi-sample VCF-file with variants to merge.')
	parser_merge.add_argument('-ploidy', metavar='PLOIDY', required=True, type=int, help='ploidy of the samples.')
	parser_merge.add_argument('-chromosomes', metavar='CHROMOSOMES', default='', help='comma separated list of chromosomes. Only output variants on these chromosomes.')
	
	parser_combine = subparsers.add_parser('combine_columns', help='combine single haplotype columns into diploid genotype columns.')
	parser_combine.add_argument('-vcf', metavar='VCF', required=True, help='VCF-file with one column per haplotype.' )
	parser_combine.add_argument('-samples', metavar='SAMPLES', help='file containing one line per sample to be merged (sample_name,h0,h1)', required=True)	
	args = parser.parse_args()

	if args.subparser_name == 'merge':
		chromosomes = None
		if args.chromosomes != '':
			chromosomes = args.chromosomes.split(',')
		run_merge(args.r, args.vcf, args.ploidy, chromosomes)
	if args.subparser_name == 'combine_columns':
		run_combine_columns(args.vcf, args.samples)
