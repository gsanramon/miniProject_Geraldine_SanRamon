import Bio
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "gsanramon@luc.edu"
import os
import gzip
import sys

option = sys.argv[1]

#define variables for SRR numbers
srr = ['SRR5660030.1', 'SRR5660033.1', 'SRR5660044.1', 'SRR5660045.1']
srr_dict = {'SRR5660030.1':'2dpi', 'SRR5660033.1':'6dpi', 'SRR5660044.1':'2dpi', 'SRR5660045.1':'6dpi'}
#to use for #4 logfile
srr_dict_2 = {'SRR5660030.1':'Donor 1 (2dpi)', 'SRR5660033.1':'Donor 1 (6dpi)', 'SRR5660044.1':'Donor 3 (2dpi)', 'SRR5660045.1':'Donor 3 (6dpi)'}

#1. Get the transcriptome files from SRA. This is skipped if running the test. Short files for are provided.

if option == "full":
  #get the files from SRA
  os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1')
  os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')
  os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')
  os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')
  
  #split the files
  os.system('fastq-dump -I --split-files SRR5660030.1')
  os.system('fastq-dump -I --split-files SRR5660033.1')
  os.system('fastq-dump -I --split-files SRR5660044.1')
  os.system('fastq-dump -I --split-files SRR5660045.1')



#2. Build the index with kallisto using the CDS features of HCMV (NCBI accession EF999921)

#retrieve the HCMV genbank file
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text")
#read the genbank file  
record = SeqIO.read(handle, "genbank")
#open output fasta file to put CDS features to be extracted from genbank format 
output_handle = open("EF999921.fasta", "w")
#initialize CDS count
cds_count = 0
#count & extract CDS as it loops through the genbank record
for feature in record.features:
  if feature.type == "CDS":
    cds_count += 1
    #Get feature qualifiers "protein_id" and "product" as identifers for the sequences to be written in fasta file
    protein_id = feature.qualifiers["protein_id"]
    product = feature.qualifiers["product"]
    feature_seq = feature.extract(record.seq)
    #FASTA output without line wrapping:
    output_handle.write(">" + ''.join(protein_id) +' '+ ''.join(product) +' from '+ str(record.id) + "\n" + str(feature_seq) + "\n")
output_handle.close()

#open logfile to write
logfile = open('miniProject.log', 'w+') 
logfile.write('The HCMV genome (EF999921) has '+ str(cds_count) +' CDS.'+'\n')

#build a transcriptome index for HCMV with kallisto
os.system('time kallisto index -i index.idx EF999921.fasta')


#3. Quantify the TPM of each CDS in each transcriptome using kallisto

#loop to quantify & create sample_info text file for sleuth (used srr list & dictionary above)
#open file to write sample info table needed for sleuth
sample_info = open('sample_info.txt', 'w+')
#write the header for sample_info
sample_info.write('sample'+'\t'+'condition'+'\t'+'path'+'\n')
#create results folder to put the 
os.mkdir('results/')
#run kallisto for each transcriptome, put into results folder & write to sample_info table
for s in srr:
  os.system('time kallisto quant -i index.idx -o results/'+s+' -b 30 -t 2 '+s+'_1.fastq '+s+'_2.fastq')
  sample_info.write(s+'\t'+srr_dict[s]+'\t'+'results/'+s+'\n')
sample_info.close()

#run sleuth.R
os.system('Rscript sleuth.R')
#os.system('nohup Rscript sleuth.R &')

#read sleuth_log.txt & write to logfile
with open("sleuth_log.txt") as lines:
    for line in lines:
        logfile.write(line)


#4. Using bowtie2 create an index for HCMV 

#build index for EF999921
os.system('bowtie2-build EF999921.fasta EF999921_index')

#loop to assemble transcriptome reads using bowtie2
for s in srr:
  os.system('bowtie2 -x EF999921_index -1 '+s+'_1.fastq -2 '+s+'_2.fastq -S '+s+'.sam --al-conc-gz '+s+'_mapped_%.fq.gz')

#write to logfile the number of reads in each transcriptome before and after the Bowtie2 mapping
for s in srr:
  pre, post = 0, 0
  for rec in SeqIO.parse(s+"_1.fastq", "fastq"):
    pre += 1
  with gzip.open(s+"_mapped_1.fq.gz","rt") as handle:
    for rec in SeqIO.parse(handle, "fastq"):
      post += 1
  #print(srr_dict_2[s]+' had '+str(pre)+' read pairs before Bowtie2 filtering and '+str(post)+' read pairs after.')
  logfile.write(srr_dict_2[s]+' had '+str(pre)+' read pairs before Bowtie2 filtering and '+str(post)+' read pairs after.'+'\n')

#5.Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. 
spades_cmd = 'spades -k 127 -t 2 --only-assembler -1 SRR5660030.1_mapped_1.fq.gz -2 SRR5660030.1_mapped_2.fq.gz -1 SRR5660033.1_mapped_1.fq.gz -2 SRR5660033.1_mapped_2.fq.gz -1 SRR5660044.1_mapped_1.fq.gz -2 SRR5660044.1_mapped_2.fq.gz -1 SRR5660045.1_mapped_1.fq.gz -2 SRR5660045.1_mapped_2.fq.gz -o SRR_assembly/'

#run SPAdes to assemble the transcriptome
os.system(spades_cmd)

#write the spades command to logfile
logfile.write(spades_cmd+'\n')

#6. Calculate the number of contigs with a length > 1000
#7. Calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length)

#go thru the input file to get contigs > 1000
#assign variable for input file after assembly via SPades
file = 'SRR_assembly/contigs.fasta'
#initialize variables
con_total_len = 0
count_contig = 0
max_conlen = 0
max_contig = str
#loop through contigs.fasta to get contig length, count contigs & total length of contigs > 1000
for record in SeqIO.parse(file, "fasta"):
    con_len = len(record)
    if con_len > 1000:
        count_contig += 1
        con_total_len += con_len
        if con_len > max_conlen:
            max_id = record.id
            max_contig = record.seq
            max_conlen = con_len
    
#write to log file the following:
#6. number of contigs with a length > 1000
logfile.write('There are '+str(count_contig)+' contigs > 1000 bp in the assembly.'+'\n')
#7. length of the assembly (the total number of bp in all of the contigs > 1000 bp in length)
logfile.write('There are '+str(con_total_len)+' bp in the assembly.'+'\n')

#write the longest contig to a file as an input for blast
contigout = open('contigout.fasta', 'w+') 
contigout.write('>'+ max_id + '\n' + str(max_contig) + '\n')
contigout.close()

#8. Use the longest contig as blast input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily

#enter the input/output files 
query_file = 'contigout.fasta'
output_file = 'blastn_results.csv'
db_name = 'local_db'
blast_log_file = 'blast_log.csv'

#make a local database only if option for args is full
if option == "full":
  #Get Betaherpesvirinae nucleotide sequences from NCBI if running full and save as Betaherpesvirinae.fasta
  input_file = 'Betaherpesvirinae.fasta'
  #make local database
  makeblast_cmd = 'makeblastdb -in '+input_file+' -out '+db_name+' -title '+db_name+' -dbtype nucl'
  os.system(makeblast_cmd)

#blastn query
blast_cmd = 'blastn -query '+query_file+' -db '+db_name+' -out '+output_file+' -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_cmd) 

#assign headers
headers = ['Subject accession', 'Percent identity', 'Alignment length', 'Start of alignment in query', 'End of alignment in query', 'Start of alignment in subject', 'End of alignment in subject', 'Bit score', 'E-value', 'Subject Title']

#write to log file the headers
logfile.write(str('\t'.join(headers))+'\n')

#open & read the blast results file generated by the blast query to write to logfile
blast_results = open(output_file, 'r')
results = blast_results.readlines()
#write to log file the top 10 hits of blast query
for line in results[0:10]:
  #print(line)
  logfile.write(str(line).replace(',','\t'))



#close logfile
logfile.close()
