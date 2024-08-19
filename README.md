The following describes code to download the data analysed, the tools required, and how and in which order to execute the scripts.

Tools
- SRA Toolkit (‘fastq-dump’ ())
- Fastqc
- Trimmomatic
- STAR
- featureCounts
- BioMart tool
- Load packages in individual scripts

Analysis steps

0 Download raw data
- Raw data in .fq format can be downloaded from SRA ( ) using the accession numbers () of the datasets.
- GTF files .....

1 Raw data preprocessing
- Check raw read sequences for quality with FASTQC version 0.11.9.
$ fastqc -o ~/path/to/outputdir/ -t 6 *.fq
- Trim off sequencing adapters, low-quality reads and short reads with Trimmomatic version 0.39 with the following parameters: ILLUMINACLIP: TruSeq2 (SE) or TruSeq3 (PE) or Nextera (PE) LEADING: 3 TRAILING: 3 SLIDINGWINDOW: 4:20 MINLEN: 20. 
$ trimmomatic PE <READ1> <READ2> -baseout Trim_Reads/<basename> ILLUMINACLIP:<FASTAadapterFILE> LEADING: 3 TRAILING: 3 SLIDINGWINDOW:4:20 MINLEN:20
- Align reads to their respective reference genomes and gene annotation files (Ensembl Protists release 51) using STAR version 2.7.8a with the following parameters: -sjdbOverhang (“read length” – 1) - limitBAMsortRAM 50000000000. 

First generate genome indices:
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir </path/to/store/genome_indices> \
--genomeFastaFiles </path/to/FASTA_file> \
--sjdbGTFfile </path/to/GTF_file> \
--sjdbOverhang <readlength -1> \
--limitBAMsortRAM 50000000000

Then align:
STAR --genomeDir </path/to/genome_indices> \
--runThreadN 6 \
--readFilesIn </path/to/FASTA_file> \
--outFileNamePrefix </path/to/store/output> \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

- Quantify the aligned reads that mapped to exonic regions with ’featureCounts’ from the Subread package version 2.0.1 with the following parameters: -s Stranded (Reverse) or Unstranded -p (for PE reads) -t exon -g gene id.


2 Merge data sets
- Transform P. falciparum and P. vivax gene ids by orthology into P. berghei gene ids using the BioMart tool in the Ensembl Protists genome browser (release 51).
- Integrate read counts for P. berghei into a single matrix with those of P. falciparum and P. vivax. 
- Retain genes with 1:1:1 orthologues across the three species.

3 Quality control
- Perform Spearman correlation test with hierarchical clustering on the 500 most variable genes to evaluate the correlation of the P. berghei datasets with the P. falciparum and P. vivax datasets. 
- Perform dimensionality reduction by principal component analysis (PCA) on the 500 most variable genes. 
- Remove batch effects in the datasets with ‘removeBatchEffect’, implemented in limma.

4 Differential gene expression analysis 

5 Coexpression clustering analysis

6 Motif activity inference

7 Cross-correlation analysis







