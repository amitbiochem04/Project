https://www.biostars.org/p/174376/
https://f1000research.com/articles/4-1080/v1
http://biocluster.ucr.edu/~rkaundal/workshops/R_feb2016/ChIPseq/ChIPseq.html

https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
http://www.bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
http://guangchuangyu.github.io/2014/01/bug-of-r-package-chippeakanno/
http://www.genomebiology.com/2015/16/1/237
http://biorxiv.org/content/early/2016/01/30/038281.1
http://biorxiv.org/content/early/2015/08/14/024620
http://biorxiv.org/content/early/2015/06/29/021642
http://nservant.github.io/HiC-Pro/MANUAL.html
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4556012/
https://github.com/ENCODE-DCC

####importnat to read 
http://genomicsclass.github.io/book/pages/classes.html

http://www.nature.com/nbt/journal/v34/n2/full/nbt.3477.html
http://www.nature.com/nbt/journal/v34/n2/full/nbt.3450.html
http://www.nature.com/nbt/journal/v34/n2/full/nbt.3468.html
Science 351, 1454–1458. doi:10.1126/science.aad9024

###TAD
http://biorxiv.org/content/early/2016/03/15/042952.full.pdf
https://www.ncbi.nlm.nih.gov/pubmed/25732821


#######analysis workflow
 
 1) Reads preprocessing
        Quality filtering (trimming)
        FASTQ quality report
 2) Alignments: Bowtie2 or rsubread
 3) Alignment stats
 4) Peak calling: MACS2, BayesPeak
 5) Peak annotation with genomic context
 6) Differential binding analysis
 7) GO term enrichment analysis
 8) Motif analysis
 
#######comand linef for the analysis 
#######Read pre-processing in R 
library(systemPipeR)
library(systemPipeRdata)
genWorkenvir(workflow="chipseq")
setwd("chipseq")
####write sample name in txt file with "FileName" "SampleName" "Factor" "SampleLong" "SampleReference" example 
            FileName SampleName Factor SampleLong SampleReference
## 1 ./data/SRR446027_1.fastq        M1A     M1  Mock.1h.A                
## 2 ./data/SRR446028_1.fastq        M1B     M1  Mock.1h.B                
## 3 ./data/SRR446029_1.fastq        A1A     A1   Avr.1h.A             M1A
## 4 ./data/SRR446030_1.fastq        A1B     A1   Avr.1h.B             M1B


targetspath <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]

#####read quality and triming 

args <- systemArgs(sysma="param/trim.param", mytargets="targets_chip.txt")
filterFct <- function(fq, cutoff=20, Nexceptions=0) {
     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
     fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
 }
preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
writeTargetsout(x=args, file="targets_chip_trim.txt", overwrite=TRUE)

########FASTQ quality report

args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip_trim.txt")
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()

#####https://scilifelab.github.io/courses/ngsintro/1604/labs/chipseq
####http://homer.salk.edu/homer/introduction/update.html


######comand line for quality check 
fasqc SRR446027_1.fastq
####triming 
java -jar $TRIMMOMATIC/trimmomatic-0.30.jar SE \
     -threads 16 \
     -phred33 \
     Sample1.fastq Sapmple1_pp.fastq \
     ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 \
     LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

####removing adptor and qulaity check for paired end result 
##https://github.com/seqan/flexbar/wiki/Manual
###flexbar -r $1 -t $2/$target -n 20 -z BZ2 -m 30 -u 0 -q 28 -a /biosw/flexbar/Adapter.fa -f sanger
flexbar -r $1 -p $2 -t $3/$target -n 20 -z GZ -m 30 -u 0 -q TAIL -qt 28 -a /biosw/flexbar/Adapter.f -qf sanger -j

##### alignment in bowtie2 for single end and paired end alignmnet 

bowtie2 -x /home/amit/genome/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome -U SRR2339650.fastq -S SRR2339650.sam
bowtie2 -x /home/amit/genome/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome -1 Reads1.fastq -2 Reads2.fastq –S DNA.sam

######peak calling by macs2
##http://cbsu.tc.cornell.edu/lab/doc/CHIPseq_workshop_20150504_lecture1.pdf
###http://www.biologie.ens.fr/~mthomas/other/chip-seq-training/

######Simply filter the bam with MAPQ (mapping quality of the reads), 5 or 10 is usually reasonable: convert sam file to bam file 
samtools view -bS  csp.sam > csp.bam
samtools sort csp.bam -o csp.sorted.bam
samtools index csp.sorted.bam
########
#samtools view -b -q 10 foo.bam > foo.filtered.bam
##or if you only want the number:
#samtools view -c -b -q 10 foo.bam
###### PhantomPeakQualTools check



#####peak call for narrow peaks:
macs2 callpeak -t IP.bam -c Input.bam -n test -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs

#####for borad regions:
macs2 callpeak -t IP.bam -c Input.bam -n test --broad -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs

######
macs2 -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01

macs2 callpeak -t SRR446029_1.fastq_trim.gz_A1.bam -c SRR446027_1.fastq_trim.gz_M1.bam -n SRR446029_1.fastq_trim.gz_A1.bam_macs2 -f BAM -g 1.2e8 -B -q 0.01 --nomodel
###if replicate is there then one can merge them. 

t  ChIP-seq treatment files; REQUIRED.
-c  Control or mock data files in either BED format or any ELAND output format specified by -format option
-n  Experiment name, which will be used to generate output file names; DEFAULT: "NA"
-f  Format of tag file, "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDMULTIPET" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE"
-g  Effective genome size (this is the size of the genome considered "usable" for peak calling). It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). It is smaller than the complete genome because many regions are excluded (telomeres, highly repeated regions...)
-B  Whether or not to save extended fragment pileup at every bp into a bedGraph file. When it's on, -w, --space and --call-subpeaks will be ignored. When --single-profile is on, only one file for the whole genome is saved
-q  The q-value (minimum FDR) cutoff to call significant regions. Default is 0.01. For broad marks, one can try 0.05 as cutoff. Q-values are calculated from p-values using Benjamini-Hochberg procedure
--nomodel Whether or not to build the shifting model. If True, MACS will not build model. by default it means shifting si

########convert bam to bed file 
######Install bedClip and bedGraphToBigWig UCSC utilities first.
###one can convert the bam to bed and the use bedtools slop to extend the reads to 3' for 200bp and then feed into bedtools coverage biostar post
bamToBed -i input.bam | slopBed -i - -g genome_file_of_chr_sizes -s -r 164 | bedToBam -i - -g genome_file_of_chr_sizes > output_extended.bam






####### peak annotation by homer 
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4938642/

##findMotifsGenome.pl test_summits.bed mm10 ER_MotifOutput/ -size 200 -mask
findMotifsGenome.pl test.bed mm10 MotifOutput/ -size 200 -mask
annotatePeaks.pl test.bed mm10 >peakannotate_homer.txt
#####bam header correction 

samtools view -H 1.bam > header
samtools reheader header 2.bam > 2.fixed.bam
samtools index 2.fixed.bam
bedtools multicov -bams 1.bam 2.fixed.bam -bed test.bed

#####Preparing ChIP-seq count table
##count one bam file:
time bedtools multicov -bams ../../data/wgEncodeSydhHistoneMcf7H3k27acUcdAlnRep1.bam -bed 1000_bedtools.bed > counts_multicov.txt






