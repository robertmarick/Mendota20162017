#!/bin/sh

#############################
# 16S_analysis.sh
# Benjamin D. Peterson

# This is the code for the 16S analysis workflow
# I used in my 5M project.
#############################

#############################
# Data transfer
#############################
#mkdir ~/5M/dataRaw/16S
#sftp bpeterson26@download.biotech.wisc.edu
#cd Bioinformatics\ Resource\ Center/McMahon/180504_BTJMK
#get ME_*

#############################
# Transfer and rename data
#############################

screen -S 5M_16S_processing

mkdir ~/5M/dataEdited/16S_processing
mkdir ~/5M/dataEdited/16S_processing/fastq_files

cd ~/5M/dataRaw/16S/
for sample in *.fastq.gz
do
  cd ~/5M/dataRaw/16S/
  newName=$(echo $sample | sed 's/_S[0-9]*_L001//g' | sed 's/_001//' | sed 's/.001//')
  echo $newName
  cp $sample ~/5M/dataEdited/16S_processing/fastq_files/$newName
  cd ~/5M/dataEdited/16S_processing/fastq_files
  gunzip $newName
done

# Make list of samples to use
cd ~/5M/dataEdited/16S_processing/fastq_files
ls *fastq \
  | sed 's/_R[1-2].fastq//' \
  | sort \
  | uniq \
  > ~/5M/dataEdited/16S_processing/16S_sample_list.txt


#############################
# Count reads
#############################
screen -S 5M_16S_readCounting
cd ~/5M/dataRaw/16S/
ancillary_info=~/5M/dataEdited/16S_processing

echo -e "seqID\tforwardReads" > $ancillary_info/seq_read_count.tsv
ls *fastq.gz | while read file16S
do
  forwardCount=$(zgrep -E -c '^@HWI-M04026|^@M04026' $file16S)
  echo -e $file16S"\t"$forwardCount
  echo -e $file16S"\t"$forwardCount >> $ancillary_info/seq_read_count.tsv
done


#############################
# Set up programs
#############################
cd ~/5M/dataEdited/16S_processing/
cp ~/programs/16S/programsFor16S.zip ./
unzip programsFor16S.zip
tar -zxvf mothur.tar.gz
tar -zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -zxvf R.tar.gz
tar -zxvf Python-2.7.13.tgz
rm -f mothur.tar.gz ncbi-blast-2.6.0+-x64-linux.tar.gz R.tar.gz Python-2.7.13.tgz

#############################
# Start mothur processing
#############################
mothur/mothur
make.file(inputdir=fastq_files, outputdir=., type=fastq)
system(mv stability.files 5M.file)
make.contigs(file=5M.file, outputdir=., processors=20)
summary.seqs()
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=600, maxhomop=8)
unique.seqs(fasta=current)
count.seqs(name=current, group=current)
align.seqs(fasta=current, reference=/home/GLBRCORG/bpeterson26/references/16S/alignment/silva.bacteria.fasta, flip=T)
summary.seqs()
# Trim the alignment. Cut those that start after 6388, and end before 25300
screen.seqs(fasta=current, count=current, summary=current, start=6388, end=25300)
summary.seqs()
# Filter out the sequences to cut out ones with a period in them
filter.seqs(fasta=current, vertical=T, trump=.)
# Remake the unique sequence list with the filtered sequences
unique.seqs(fasta=current, count=current)
# Files after unique:
# /mnt/bigdata/linuxhome/bpeterson26/5M/dataEdited/16S_processing/5M.trim.contigs.good.unique.good.filter.count_table
# /mnt/bigdata/linuxhome/bpeterson26/5M/dataEdited/16S_processing/5M.trim.contigs.good.unique.good.filter.unique.fasta
# Precluster the data to denoise it. Speeds up computation
pre.cluster(fasta=current, count=current, diffs=4)
summary.seqs()
# Remove chimeras
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
quit()

#############################
# Mothur's little helper
#############################
cp 5M.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta 5M_raw.fasta
cp 5M.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table 5M_raw.count
mkdir processingFiles
mv 5M.* processingFiles
rm -r -f fastq_files
sed 's/-//g' 5M_raw.fasta > 5M_clean.fasta
sed 's/-//g' 5M_raw.count > 5M_clean.count
mv 5M_raw.fasta processingFiles
mv 5M_raw.count processingFiles

# 1546595 (1.5 million) sequences in this list
# How many are singletons?
awk -F '\t' ' $2 != 1 { print $0 }' 5M_clean.count | wc -l
# Shit, there's a lot. Let's remove those.
awk -F '\t' ' $2 != 1 { print $0 }' 5M_clean.count > 5M_clean_trimmed.count
# Let's download this file and check it out in R.

#############################
#############################
# Run TaxAss
#############################
#############################

cp ~/programs/16S/taxAssScripts.zip ./
unzip taxAssScripts.zip
chmod +x *.R *.py

#############################
# Make blast DB
#############################
# This only needs to be run once
references=~/references/16S/
ncbi-blast-2.6.0+/bin/makeblastdb -dbtype nucl \
                                  -in $references/FW.fasta \
                                  -input_type fasta \
                                  -parse_seqids \
                                  -out $references/FW.fasta.db

#############################
# Run BLAST on samples
#############################
# Blast OTUs against FW blast database
cd ~/5M/dataEdited/16S_processing/
references=~/references/16S/
ncbi-blast-2.6.0+/bin/blastn -query 5M_clean.fasta \
                            -task megablast \
                            -db $references/FW.fasta.db \
                            -out OTU.custom.blast \
                            -outfmt 11 \
                            -max_target_seqs 5
# Reformat the blast results
ncbi-blast-2.6.0+/bin/blast_formatter -archive OTU.custom.blast \
                                      -outfmt "6 qseqid pident length qlen qstart qend" \
                                      -out OTU.custom.blast.table
# Use Robin's script to calculate the full length pident of each blast result.
Rscript calc_full_length_pident.R OTU.custom.blast.table OTU.custom.blast.table.modified

#############################
# Pull out sequences to run with Silva vs. FreshTrain
#############################
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.above.98 98 TRUE
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.below.98 98 FALSE
python find_seqIDs_blast_removed.py 5M_clean.fasta OTU.custom.blast.table.modified ids.missing

cat ids.below.98 ids.missing > ids.below.98.all
python create_fastas_given_seqIDs.py ids.above.98 5M_clean.fasta otus.above.98.fasta
python create_fastas_given_seqIDs.py ids.below.98.all 5M_clean.fasta otus.below.98.fasta

#############################
# Classify sequences
#############################
mothur/mothur "#classify.seqs(fasta=otus.above.98.fasta, template=/home/GLBRCORG/bpeterson26/references/16S/FW.fasta, taxonomy=/home/GLBRCORG/bpeterson26/references/16S/FW.taxonomy, method=wang, probs=T, processors=18)"
mothur/mothur "#classify.seqs(fasta=otus.below.98.fasta, template=/home/GLBRCORG/bpeterson26/references/16S/silva.fasta, taxonomy=/home/GLBRCORG/bpeterson26/references/16S/silva.taxonomy, method=wang, probs=T, processors=18)"
cat otus.above.98.FW.wang.taxonomy otus.below.98.silva.wang.taxonomy > 5M_clean.taxonomy
sed "s/[[:blank:]]/,/" 5M_clean.taxonomy | \
  sed "s/;/,/g" | \
  sed 's/([0-9]*)//g' | \
  sed 's/,$//g' \
  > 5M_cleaner.taxonomy



#############################
# Calculate % abundance of 1 read in the smallest metagenome
#############################

# Count up reads in each sample

cd ~/5M/dataEdited/16S_processing/fastq_files
for sample in *fastq
do

  echo "counting" $sample
  lineCount=`wc -l < $sample`
  echo -e $sample"\t"$lineCount >> fastq_line_counts_not_reads.tsv

done

# Download to my computer,




#############################
# Clean up!
#############################
rm -r mothur ncbi-blast-2.6.0+/ programsFor16S.zip R Python-2.7.13/ __MACOSX/ taxAssScripts.zip
mv mothur*logfile OTU* otus* *.R *.py ids* processingFiles


# How many sequences were removed at each step?
cd ~/HellsCanyon/dataEdited/16S_processing/processingFiles
wc -l *fasta

#613570 HCC_raw.fasta
#      0 HCC.scrap.contigs.fasta
#5819806 HCC.trim.contigs.fasta
#4743302 HCC.trim.contigs.good.fasta
#3457514 HCC.trim.contigs.good.unique.fasta
#3371748 HCC.trim.contigs.good.unique.good.filter.fasta
#2042824 HCC.trim.contigs.good.unique.good.filter.unique.fasta
# 695634 HCC.trim.contigs.good.unique.good.filter.unique.precluster.fasta
# 613570 HCC.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
#21357968 total

# Looks like sequence list was trimmed the most by the preclustering.
