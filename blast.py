#! /bin/bash/env python
from Applications import NcbiblastpCommandline
import math
import os
from decimal import *
import settings
added = []
class Blaster(object):
	def __init__(self):
		pass
		
	def blast(self, evalue):
		# create a database and conduct the all-vs-all BLAST search against it
		print 'Creating database containing all sequences from all genomes...'
		os.system('makeblastdb -in venninator_combined.fa -dbtype prot -out venninator_database')
		combinedfile_splits = [splitfile for splitfile in os.listdir('.') if splitfile.startswith('venninator_combined_')]
		pctdone = 0
		if settings.allvsallfile == '':
			print 'Executing all-vs-all BLAST search (this will take a while) ...'
			for splitfile in combinedfile_splits:
				blastsearch = NcbiblastpCommandline(cmd='blastp', query=splitfile, db='venninator_database', outfmt='"6 qseqid sseqid qcovs evalue"', evalue=settings.evalues[0], max_hsps_per_subject='1', out='out_blastp_' + settings.evalues[0] +'_' + splitfile.split('_')[2] + '.tmp', num_threads=settings.numthreads)
				blastsearch()
				pctdone += 10
				print '%d%% done ...' % pctdone
		else:
			with open('out_blastp_' + settings.evalues[0] + '_combined.tmp', 'w') as copyfile:
				with open(settings.allvsallfile, 'r') as avafile:
					for line in avafile:
						line = line.strip()
						copyfile.write(line + '\n')

	def trimmer(self, evalue):
		# extract into new file "query, top hit, e-value" from rows that meet conditions (no self-hits; consolidate reciprocal hits; >70% coverage).
		print 'Taking all-vs-all BLAST output and trimming (removing duplicates, length below threshold, combining e-values)...'
		blastoutput  = open('../out_blastp_' + settings.evalues[0] + '_combined.tmp')
		outfile = open('out_blastp_' + evalue + '_trimmed.tmp', 'w')
		seen = set()
		withevalue = {}
		# This generator expression goes line by line ONCE (instead of iteratively),
		# and returns the values of each line, excluding self-hits (A1 -> A1) and hits with length less than 70%.
		line = (line.strip().split('\t')[0:4] for line in blastoutput if len(line) > 0 and line.split('\t')[0] != line.split('\t')[1] and int(line.split('\t')[2]) >= settings.length and Decimal(line.split('\t')[3]) < Decimal(evalue))
		for i in line:
			pair = ','.join(sorted(i[0:2]))		# Use sort() in order to make sure values of A1 -> B2 and B2 -> A1 are combined.
			if pair not in seen:
				seen.add(pair)					# Seen is the master set of pairs
				withevalue[pair] = [i[3]]		# if an A->B not in seen, then we define that A->B as the key, and create a list with this first e-value for its value
			elif pair in seen:					# If we've already added this A->B pair,
				withevalue[pair].append(i[3])	# then we append the e-value (value) to this existing pair (key) in our dictionary.
		for x in withevalue:
			outfile.write(x + ',' + ','.join(withevalue[x]) + '\n')
		outfile.close()
	
	def homolog_splitter(self, evalue):
		# separate homologs among species from in-paralogs within species
		# this is done such that their weights can be normalized respectively
		print 'Separating orthologs and paralogs...'
		with open('out_blastp_' + evalue + '_trimmed.tmp', 'r') as infile:
			for line in infile:
				line = line.strip().split(',')
				qid = line[0].split('_',1)[0]
				hid = line[1].split('_',1)[0]
				line = ','.join(line)
				with open('out_blastp_' + evalue + '_trimmed_splitted.tmp', 'a') as outfile:
					if qid == hid:
						outfile.write('p,' + line + '\n')
					elif qid != hid:
						outfile.write('o,' + line + '\n')
			
	def recip_averager(self, evalue):
		outfile = open('output_readyformcl_' + evalue + '.tmp', 'w')
		with open('out_blastp_' + evalue + '_trimmed_splitted.tmp', 'r') as infile:
			finalaveraged = {}
			for line in infile:
				line = line.strip().split(',')
				pair = '\t'.join(line[1:3])
				finalaveraged[pair] = ''
				justevalues = line[3:]
				eneglogs = []
				for evalue in justevalues:
					if evalue == '0' or evalue == '0.0' or evalue == 0 or evalue == 0.0:
						evalue = 308.0	# Rounding up -log10(2.225074e-308), the lowest possible e-value before 0.0
					else:	
						evalue = -math.log10(Decimal(evalue))
					eneglogs.append(evalue)
				averaged = sum(eneglogs) / len(eneglogs)
				finalaveraged[pair] = str(averaged)
			for x in finalaveraged:
				outfile.write(x + '\t' + finalaveraged[x] + '\n')
