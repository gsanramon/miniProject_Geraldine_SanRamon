# miniProject_Geraldine_SanRamon
Documents included in this repo are the following:
1. miniTest_copy folder that has the following files:

  a. run_me.py - the python script to run with the option to run the full or short test version.
  
      So in the command line, please type the following if running the full code:
      $ python3 run_me.py full
      
      If running the test option write the following in the command line:
      $ python3 run_me.py test
      
     Note: If running the full run, you will need to get the nucleotide sequences of Betaherpesvirinae subfamily and save it as Betaherpesvirinae.fasta in the same working directory. This Betaherpesvirinae.fasta is ~140Mb so it's not included here. 
     
  b. sleuth.R - Rscript that will be called in the python script
  
  c. local_db files (3 total) - this is the local blast database already generated if running the test option only.
  
  d. fastq files (8 total) - first 10,000 of paired-end reads of th 4 transcriptomes downloaded from SRA. The full reads will be downloaded as part of the full run.
  
2. logfile - logfile after my full run of the python script
3. README.md - this file you're reading now

The following software must be installed to run the code:
1. kallisto - used to build index & quantify TPM
2. R - R package sleuth must also be installed in addition to R
3. Bowtie2 - mapping software used to create index to use for assembly via SPAdes
4. SPAdes - de novo genome assembly
