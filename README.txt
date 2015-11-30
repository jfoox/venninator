// VENNINATOR v0.1 (November 30, 2015)
// Gene Content Matrix Assembler using Orthologs Clustered under Expect Value Gradients
// Jonathan Foox et al. - http://www.github.com/jfoox/venninator
// See publication in XXX
// License: MIT License (MIT)

==========================
Download and installation:
==========================
Venninator is extremely simple to download and use. All required Python scripts are 
included in the Venninator package. No compilation is required. Note that Venninator is
intended to be used with a Linux/Mac OS X enviroment, and is not supported for Windows.

A pre-compiled version of MCL is also included, version mcl-14-137. You may manually
download and include a newer version of MCL if desired (http://micans.org/mcl/). Assuming
that the command line usage of MCL is not changed, then nothing needs to be changed within
the Venninator pipeline. (If change is required, edit line 126 of venninator.py.)

=============
Requirements:
=============
* Folder containing all FASTA files of interest
* FASTA files must end with extension '.fa'
* Venninator only processes protein-coding peptide sequences (no nucleotide sequences)

* The longest step of the process is the initial BLAST search, in which every sequence
  from every FASTA file is compared against every other sequence (all-vs-all). If you
  would rather conduct this search ahead of time (for instance on a high-powered cluster), 
  you must do the following:
    * Ensure that all sequences comprise the BLAST database, and that all sequences are
      used as the query
    * Conduct BLASTP with the following options:
      - outfmt='6 qseqid sseqid qcovs evalue'
      - evalue=LOWEST_EVALUE_YOU_WISH_TO_USE (e.g., 1e-5)
      - max_hsps_per_subject=1
    * Use the --allvsall (-z) flag when loading Venninator to specify the location of the
      output file
    * Ensure that, when specifying evalues using --evalues, the lowest value = the value
      you used for the BLASTP search

=============
Dependencies:
=============
* Python 2.7
* A local installation of BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* Included in package is the required library of BioPython (http://biopython.org/wiki/Main_Page)
* Included in package is pre-compiled latest version of MCL

===========
How to run:
===========
To see help for all available options:
[user folder]$ python venninator.py -h

Example run:
[user folder]$ python venninator.py --annotate --cleanup --directory path/to/folder --evalues 5 25 75 125 200 --length 80 --numthreads 16

