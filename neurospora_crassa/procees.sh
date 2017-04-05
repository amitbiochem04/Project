
#####Map against the reference genome using bowtie2
http://www.ngscourse.org/Course_Materials/rna_seq/tutorial_htseq/rna_seq_htseq.html
bowtie2-build Neurospora_crassa.NC12.dna.toplevel.fa neurospora_crassa_bowtie_index
tophat -p 4 -G /Users/amitsingh/Desktop/genome/neurosopra/Neurospora_crassa.NC12.34.gtf -o test1 --transcriptome-index= /Users/amitsingh/Desktop/genome/neurosopra/neurospora_crassa_bowtie_index lane5Dfrq0_sequence.fastq

htseq-count -f bam -m intersection-nonempty -r pos --type=exon --idattr=gene_id --stranded=yes accepted_hits.bam /Users/amitsingh/Desktop/genome/neurosopra/Neurospora_crassa.NC12.34.gtf >test.txt

###single end 
featureCounts -a /Users/amitsingh/Desktop/genome/neurosopra/Neurospora_crassa.NC12.34.gtf -t exon -g gene id -o counts.txt mapping accepted_hits.bam

featureCounts -Q 10 -M -s 0 -T 1 -a Homo_sapiens.GRCh38.76.gtf -o test.bam.featureCount.cnt test.bam
featureCounts -Q 10 -M -s 0 -T 8 -a Homo_sapiens.GRCh38.76.gtf -o test.bam.featureCount.cnt test.bam

http://staff.um.edu.mt/jebej02/blog/2016/05/23/installing-cufflinks-rnaseq-on-ubuntu/
