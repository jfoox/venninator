#! /bin/bash/env python
# encoding: utf-8
# please refer to README.md for instructions

# import dependencies
import settings
from seqcurator import *
from blast import Blaster
from postmcl import PostMCL
import argparse
import time
import os
import shutil
import sys
settings.init()

header = '''\
**************************************************
*********** Welcome to VENNINATOR v0.1 ***********
**************************************************
   Gene Content Matrix Assembler using Orthologs
     Clustered under Expect Value Gradients.'''

### command line input
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=(header))
parser.add_argument('-a', '--annotate', action='store_true', help='Generate annotation file for each gene content matrix at each e-value specified')
parser.add_argument('-c', '--clean', action='store_true', help='Full clean-up, leaving only NEXUS and annotation files')
parser.add_argument('-d', '--directory', metavar='', default='fasta_folder', help='Folder containing your FASTA files, if you have renamed it (default: fasta_folder)')
parser.add_argument('-e', '--evalues', metavar='', nargs='+', type=int, default=[5,25,75], help='E-values to use for analyses; write as integers separated by  (default: 5 25 75)')
parser.add_argument('-l', '--length', metavar='', type=int, default=70, help='Minimum percentage of sequence\'s length recovered that is required to retain BLAST match (default: 70)')
parser.add_argument('-t', '--numthreads', metavar='', type=int, default=8, help='Number of threads to use for BLAST searches (default: 8)')
parser.add_argument('-z', '--allvsallfile', metavar='', default='', help='Location of all-vs-all BLAST output file if done separately (see documentation)')
args = parser.parse_args()
print header

### interpret user's command line input and ensure they want to proceed
# all-vs-all output
if args.allvsallfile != '':
	if not os.path.isfile(args.allvsallfile) and not os.path.isfile(os.getcwd() + '/' + args.allvsallfile):
		sys.exit('%s not a valid file. Please try again.' % args.allvsallfile)
	else:
		settings.allvsallfile = os.path.abspath(args.allvsallfile)

# working directory
if args.directory != '':
	if not os.path.exists(args.directory):
		sys.exit('%s not a valid directory. Please make sure you specify a folder within the parent Venninator folder.' % args.directory)
	else:
		os.chdir(args.directory)

# e-value list
curator = Curate()
settings.user_evalues = args.evalues
curator.evalue_getter()
first_evalue = settings.evalues[0]

# remainder
settings.length = args.length
settings.numthreads = args.numthreads

### shall we proceed?
print '''
--------------------------------------------------
 Venninator will run with the following settings:
--------------------------------------------------'''
print 'Working directory: %s' % os.getcwd()
curator.filesgetter()
print 'E-values: %s' % ', '.join(settings.evalues)
if args.allvsallfile == '':
	print 'Minimum percentage of sequence\'s length recovered that is required to retain BLAST match: %s' % args.length
	print 'Number of threads to use for all-vs-all BLAST search: %s' % args.numthreads
else:
	print 'All-vs-all output file: %s' % settings.allvsallfile
if args.clean:
	print 'Full cleanup: Yes'
else:
	print 'Full cleanup: No'
if args.annotate:
	print 'Annotate matrices: Yes'
else:
	print 'Annotate matrices: No'
while True:
	decision = raw_input('If the above settings look correct, please press Y to proceed. Otherwise, press Q to quit. ')
	if decision == 'q' or decision == 'Q':
		sys.exit()
	elif decision == 'y' or decision == 'Y':
		break
	else:
		print 'Invalid selection, please try again.'
		continue
		
##########################################################################################
######################################## ANALYSIS ########################################
##########################################################################################
### file curation 
print 'Creating unique identifiers for each genome file ...'
curator.unique_id()
print 'Concatenating genome files ...'
curator.concatenation('venninator_combined.fa', settings.filelist, '_uniqueid.tmp', True)
print 'Splitting concatenated genomes file for parallel processing ...'
split = Splitter()
split.splitter('venninator_combined.fa', 10)
print 'Proceeding with BLAST analysis ...'
blaster = Blaster()
blaster.blast(first_evalue)
if args.allvsallfile == '':
	print 'Sorting all-vs-all BLAST output ...'
	settings.outblastp_splitfiles = [outblastpfile for outblastpfile in os.listdir('.') if outblastpfile.startswith('out_blastp_' + first_evalue + '_a')]
	curator.concatenation('out_blastp_' + first_evalue + '_combined.tmp', settings.outblastp_splitfiles, '', False)

### analysis for each e-value
for evalue in settings.evalues:
	print 'Now working at e-value %s ...' % evalue
	evalue_folder = 'output_' + evalue
	if not os.path.exists(evalue_folder):
		os.makedirs(evalue_folder)
	os.chdir(evalue_folder)

	print 'Pruning all-vs-all output at this e-value ...'
	blaster.trimmer(evalue)
	print 'Identifying out-orthologs and in-paralogs ...'
	blaster.homolog_splitter(evalue)
	print 'Calculating final similarity value for each pairwise comparison ...'
	blaster.recip_averager(evalue)
	print 'Generating MCL clusters ...'
	os.system('mcl output_readyformcl_' + evalue + '.tmp --abc -I 1.5 -o venninator_clusters_' + evalue + '.mcl')
# ../.././mcl
	print 'Generating binary matrix in NEXUS format ...'
	postmcl = PostMCL()
	postmcl.matrixmaker(evalue)
	print 'Annotating each cluster ...'
	if args.annotate:
		postmcl.annotater(evalue)
	
	if args.clean:
		print 'Cleaning subfolder ...'
		tmpremove = [tmp for tmp in os.listdir('.') if tmp.endswith('.tmp') or tmp.endswith('.mcl') or tmp.startswith('annotation_tobeblasted_')]
		for tmp in tmpremove:
			os.remove(tmp)
	os.chdir('..')

### final cleanup
if args.clean:
	for i in range(0, (settings.flatfiles_foldernumber + 1)):
		shutil.rmtree('flatfiles' + str(i))
	restremove = [rest for rest in os.listdir('.') if rest.endswith('.tmp') or rest.startswith('venninator_') or rest.startswith('swissprot')]
	for rest in restremove:
		os.remove(rest)
print 'Done!'
