#! /bin/bash/env python
from Applications import NcbiblastpCommandline
from seqcurator import *
import gzip
import math
import random
import re
import urllib

class PostMCL(object):
	def __init__(self):
		pass
		
	def matrixmaker(self, evalue):
		header = False
		for key in settings.species_id.keys():
			key = key.strip()
#			value = settings.species_id.values()
			matrixline = []
			with open('venninator_clusters_' + evalue + '.mcl', 'r') as clusterfile:
				for line in clusterfile:
					if key in line:
						matrixline.append('1')
					elif key not in line:
						matrixline.append('0')
			matrixline = ''.join(matrixline)
			with open('matrix_' + evalue + '.nex', 'a') as matrixfile:
				if header == False:
					matrixfile.write('#NEXUS\nBegin data;\nDimensions ntax=' + str(len(settings.species_id)) + ' nchar=' + str(len(matrixline)) + ';\nformat symbols=\"01\";\nMatrix\n')
					header = True
				matrixfile.write(settings.species_id[key] + '\t' + matrixline + '\n')
		with open('matrix_' + evalue + '.nex', 'a') as matrixfile:
			matrixfile.write(';\nEnd;')
				
	def annotater(self, evalue):
		linenum = 0
		### ANNOTATER STEP ONE: 
		#### Create master fasta file as "Genename\tSEQUENCE".
		### Like creating deinterleaved file, but with tabs instead of linebreaks.
		### (This way, when we find the gene name, we don't have to use a memory-intensive index step
		### to retrieve the sequence on the next line; instead, we split by tab, and output each element
		### as its own line in the outfile.)
		# Only do this one time during first e-value...
		if not os.path.isfile('../venninator_combined_forannotation.tmp'):
			with open('../venninator_combined_forannotation.tmp', 'w') as outfile:
				for name in settings.filelist:
					with open('../' + name + '_uniqueid.tmp') as infile:
						for line in infile:
							line = line.strip()
							if len(line) > 0:
								if line[0] == '>':
									line = line[1:]						# Removing the '>', which we will add to the final file; this makes equivalency checks below much easier
									if linenum == 0:					# So the first line isn't just a carriage return
										outfile.write(line + '\t')
										linenum += 1
									else:
										outfile.write('\n' + line + '\t')
								else:
									outfile.write(line)
						outfile.write('\n')						   # adds a linebreak at end of file such that 1st line of next file doesn't overlap with last line of previous

		### ANNOTATOR STEP TWO: 
		### From master file above, create sub-fasta file of sequences to be searched against SWISS-PROT database.
		print 'Retrieving sequences from clusters to be annotated ...'
	
		# Populate a dictionary with all sequences from master file
		genes_all = {}
		with open('../venninator_combined_forannotation.tmp', 'r') as genes_all_file:
			for line in genes_all_file:
				if '\t' in line:
					line = line.strip().split('\t')
					genes_all[line[0]] = line[1]
				
		# Generate the genes of interest you want to extract from genes_all dictionary above
		# Also populate clusters_nameandnum, so we maintain association btwn cluster numbers and genes
		genes_of_interest = set()
		clusters_nameandnum = {}
		clusternum = 1			
		with open('venninator_clusters_' + evalue + '.mcl', 'r') as genes_of_interest_file:
			for line in genes_of_interest_file:
				if '\t' in line:
					seqname = line.strip().split('\t')[0]
				else:
					seqname = line.strip()
				genes_of_interest.add(seqname)
				clusters_nameandnum[seqname] = clusternum
				clusternum += 1
		# Find the overlap between the total sequences and the desired ones
		genes_all_justnames = set(genes_all.keys())
		toextract = genes_all_justnames.intersection(genes_of_interest)
		
		# Use 'toextract' set to generate desired file
		with open('annotation_tobeblasted.tmp', 'w') as extractfile:
		    for name in toextract:
		    	extractfile.write('>' + name + '\n' + genes_all[name] + '\n')


		### ANNOTATER STEP THREE:
		### Execute BLASTP search of generated fasta file against Swiss-PROT database.
		# First, download the latest Swiss-PROT database and decompress it:
		# Do this only once ...
		if not os.path.isfile('../swissprot.tmp'):
			print 'Downloading and decompressing SWISS-PROT database ...'
			urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz', '../swissprot.fa.gz.tmp')
			with gzip.open('../swissprot.fa.gz.tmp', 'rb') as swissprot_compressed:
				with open('../swissprot.tmp', 'w') as swissprot_decompressed:
					for line in swissprot_compressed:
						swissprot_decompressed.write(line)
			swissprot_compressed.close()
		split = Splitter()
		split.splitter('annotation_tobeblasted.tmp', 10)
	
		# Set up BLAST search database:
		if not os.path.isfile('../swissprot_database.phr'):
			os.system('makeblastdb -in ../swissprot.tmp -dbtype prot -out ../swissprot_database')

		# Do segmented BLASTs of files to be annotated versus SWISS-PROT:
		annot_blast_files = [annot for annot in os.listdir('.') if annot.startswith("annotation_tobeblasted_")]
		totalnum = len(annot_blast_files)
		currentnum = 0
		print 'Now searching sequences of each cluster against the SWISS-PROT database ...'
		for file in annot_blast_files:
			currentnum += 1
			suffix = file.split('_',-1)[2]
			outname = 'annot_blastp_output_' + suffix + '.tmp'
			swissprotsearch = NcbiblastpCommandline(cmd='blastp', query=file, db='../swissprot_database', out=outname, outfmt='"6 qseqid sseqid qcovs evalue"', evalue='1e-5', max_target_seqs='1', max_hsps_per_subject='1', num_threads=settings.numthreads)
			print 'Now searching %s of %s cluster sequence files ...' % (currentnum, totalnum)
			swissprotsearch()
		
		# Concatenate the segmented BLASTP output files
		annot_blast_outfiles = [annotout for annotout in os.listdir('.') if annotout.startswith("annot_blastp_output_")]
		with open('annot_blastp_output_all.tmp', 'w') as outfile:
			for file in annot_blast_outfiles:
				with open(file) as annotfile:
					for line in annotfile:
						line = line.strip()
						outfile.write(line + '\n')

		# Unique the concatenated output file
		concatoutput_seen = set()
		with open('annot_blastp_output_all_uniqued.tmp', 'w') as concatuniqued:
			with open('annot_blastp_output_all.tmp', 'r') as concatoutput:
				for line in concatoutput:
					line = line.strip().split('\t')
					if line[0] not in concatoutput_seen:
						concatoutput_seen.add(line[0])
						concatuniqued.write('\t'.join(line) + '\n')
		
		
		### ANNOTATER STEP FOUR: Retrieve relevant SWISS-PROT flat files.
		# Do this using several folders so that none holds > 10,000		
		# Alongside, create dictionary with folder number for each ID,
		# so it can be retrieved for annotation file.
		print 'Retrieving Swiss-PROT flat files ...'

		# Create new folder and change working directory to that folder
		def foldermaker(foldernumber):
			if not os.path.exists('../flatfiles' + str(foldernumber)):
				os.makedirs('../flatfiles' + str(foldernumber))

		with open('annot_blastp_output_all_uniqued.tmp', 'r') as full_annotfile:
			foldermaker(settings.flatfiles_foldernumber)
			for line in full_annotfile:
				line = line.strip().split('\t')
				flatfile_id = line[1].split('|',-1)[1]
				if flatfile_id not in settings.flatfiles_and_foldernums:
					settings.flatfiles_added += 1
					settings.flatfiles_and_foldernums[flatfile_id] = settings.flatfiles_foldernumber
					working_foldernumber = str(settings.flatfiles_foldernumber)
					if not os.path.exists('../flatfiles' + working_foldernumber + '/%s.uni' % flatfile_id):
						try:
							urllib.urlretrieve('http://www.uniprot.org/uniprot/' + flatfile_id + '.txt', '../flatfiles' + working_foldernumber + '/%s.uni' % flatfile_id)
						except Exception:
							continue	
				if settings.flatfiles_added > 9999: 
					settings.flatfiles_foldernumber += 1
					foldermaker(settings.flatfiles_foldernumber)
					settings.flatfiles_added = 0
		

		### ANNOTATER STEP FIVE: Generate the annotation file.
		# Create dictionary: 1:id1_PAC:15701388, 2:id1_PAC:15703879, etc.
		# For each line of output_all, extract relevant information (GO, KEGG, PFAM, etc.) for each cluster rep
		# When writing that line, retrieve associated key (cluster#) for that value (sequencename)
		# Once all lines are written, sort by clusternum.
		print 'Generating annotation file ...'
		
		# get all the variables from the specific UniProt flatfile.		
		def flatfile_infograbber(current_flatfile):
			flatfileID = 'No_ID'
			description = 'No_Description'
			KEGG = 'No_KEGG'
			KEGG_KO = 'No_KEGGKO'
			GO = 'No_GO'
			Pfam = 'No_Pfam'

			for line in open('../flatfiles' + str(current_flatfile_foldernum) + '/' + current_flatfile + '.uni', 'r'):
				line = line.strip()
				if line[0:3] == 'ID ':
					flatfileID = line.split(' ')[3]
				if line[0:2] == 'DE' and description == 'No_Description':
					description = line.split('=')[1][:-1]
				if line[0:2] == 'DR':
					if line[5:9] == 'KEGG':
						try:
							KEGG = line.replace(' ','').split(';')[1]
						except Exception:
							continue
					if line[5:7] == 'KO':
						try:
							KEGG_KO = line.replace(' ','').split(';')[1]
						except Exception:
							continue
					if line[5:7] == 'GO':
						line = line.strip().split(';')
						if GO == 'No_GO':
							GO = line[1]
						else:
							GO += '; %s' % line[1]
					if line[5:9] == 'Pfam':
						Pfam = line.replace(' ','').split(';')[1]
				
			return '%s\t%s\t%s\t%s\t%s\t%s' % (flatfileID, description, KEGG, KEGG_KO, GO, Pfam)

		# Cluster by cluster, identify the corresponding number, name,
		# and then use above function to pull corresponding annotational info
		masternum = 1
		clusters_covered = []
		with open(evalue + '_annotations_unsorted.tmp', 'w') as finalannotfile:
			with open('annot_blastp_output_all_uniqued.tmp', 'r') as finalconcat:
				finalannotfile.write('Cluster #\tRepresentative Sequence\tUni-PROT ID\tDescription\tKEGG\tKEGG KO\tGO\tPfam\n')
				for line in finalconcat:
					line = line.strip().split('\t')
					current_seq = line[0]
					current_flatfile = line[1].split('|',-1)[1]
					current_flatfile_foldernum = settings.flatfiles_and_foldernums[current_flatfile]
					current_clusternum = clusters_nameandnum[current_seq]
					clusters_covered.append(current_clusternum)
					
					# Find if identified flatfile exists; if so, retrieve information
					if os.path.exists('../flatfiles' + str(current_flatfile_foldernum) + '/' + current_flatfile + '.uni'):
						restofrow = flatfile_infograbber(current_flatfile)
					else:
						restofrow = 'NULL (no annotation found)'
					
					# Write out line	
					finalannotfile.write(str(current_clusternum) + '\t' + current_seq.replace('@@', ',').replace('$$', ' ') + '\t' + restofrow + '\n')
			num_clusters = sum(1 for line in open('annot_blastp_output_all_uniqued.tmp'))
			for i in range(1, num_clusters):
				if i not in clusters_covered:
					finalannotfile.write(str(i) + '\t' + 'NULL (no annotation found)' + '\n')
					
		# Sort final annotation file by cluster number and return to results folder
		os.system('sort -n %s_annotations_unsorted.tmp > annotations_%s.txt' % (evalue, evalue))
