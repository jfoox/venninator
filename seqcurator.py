#! /bin/bash/env python
import re
import os
import math
import settings

class Curate(object):
	def __init__(self):
		pass
			
	def filesgetter(self):
		settings.filelist = [file.split('.')[0] for file in os.listdir(os.getcwd()) if file.endswith('.fa')]
		print 'FASTA files that will be analyzed: '
		filenum = 1
		for file in settings.filelist:
			print str(filenum) + '. ' + file
			filenum += 1
	
	def unique_id(self):
		# assign unique species prefix to all sequences for that species
		id = 1
		for file in settings.filelist:
			with open(file + '_uniqueid.tmp', 'w') as outfile:
				with open(file + '.fa') as infile:
					for line in infile:
						if line[0] == '>':
							line = line.replace('>', '>id%d_' % id).replace(' ', '$$').replace(',', '@@')
						outfile.write(line)
				addid = 'id' + str(id)
				settings.species_id[addid] = file
				id += 1

	def concatenation(self, concatoutput, listoffiles, suffix, fasta):
		# concatenates files stored in list
		linenum = 0
		with open(concatoutput, 'w') as outfile:
			for name in listoffiles:
				with open(name + suffix) as infile:	# The following deinterleaves fasta, so sequences can be queried easily downstream
					for line in infile:
						line = line.strip()
						if fasta == True:
							if len(line) > 0:
								if line[0] == '>':
									if linenum == 0:					# So the first line isn't just a carriage return
										outfile.write(line + '\n')
										linenum += 1
									else:
										outfile.write('\n' + line + '\n')
								else:
									outfile.write(line)
						elif fasta == False:
							if len(line) > 0:
								outfile.write(line + '\n')
					if fasta == True:
						outfile.write('\n')	   # adds a linebreak at end of file such that 1st line of next file doesn't overlap with last line of previous
				linenum = 0

	def evalue_getter(self):
		if settings.user_evalues != '':
			settings.user_evalues = ','.join(map(str, settings.user_evalues))
			settings.user_evalues = settings.user_evalues.replace(' ',',').split(',')
			settings.user_evalues = sorted(settings.user_evalues, key=int)
			for user_evalue in settings.user_evalues:
				settings.evalues.append('1e-' + user_evalue)
		else:
			settings.evalues = ['1e-5', '1e-25', '1e-75']


class Splitter(object):
	def __init__(self):
		pass
						
	def splitter(self, file, partitions):
		# divides file into subfiles for parallelization
		file_noext = file.split('.',-1)[0]
		file_justext = file.split('.',-1)[1]
		with open(file, 'r') as concatfile:
			num_lines = sum(1 for line in concatfile)
			per_split = math.ceil(((num_lines / partitions) / 2) * 2)
			if '.' in str(per_split):
				per_split = str(per_split).split('.',-1)[0]
			if per_split == 0:
				per_split = 1
			os.system('split -l %s %s.%s %s_' % (per_split, file_noext, file_justext, file_noext))
