#! /bin/bash/env python
def init():
	global user_evalues
	global allvsallfile
	global similarity
	global numthreads
	global filelist
	global species_id
	global evalues
	global outblastp_splitfiles
	global flatfiles_and_foldernums
	global flatfiles_foldernumber
	global flatfiles_added
	global current_numclusters
	user_evalues = ''
	allvsallfile = ''
	similarity = 70
	numthreads = 8
	species_id = {}
	filelist = []	
	evalues = []
	outblastp_splitfiles = []
	flatfiles_and_foldernums = {}
	flatfiles_foldernumber = 0
	flatfiles_added = 0
	current_numclusters = 0
