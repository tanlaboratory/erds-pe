#ERDS-pe
#Introduction
ERDS-pe is a tool designed to detect detect CNVs from whole-exome sequencing (WES) data. ERDS-pe employs RPKM and principal component analysis to normalize WES data and incorporates RD and single-nucleotide variation information together as a hybrid signal into a paired hidden Markov model to infer CNVs from WES data.
#Installation
ERDS-pe is easy to run. you just need:

1. Python 2.7+ (Python 3 is not yet supported)

2. Numpy, PyVCF, Pysam libraries installed.

3. Target(Probe) files in bed format.


#Running 

usage: python erds_pe.py <COMMAND> [OPTIONS] 

1. Getting RD and transformming to RPKM format

	This command is to extract the read depth (RD) signal from the BAM files and to calculate RPKM values for all samples. 

		usage: python erds_pe.py rpkm [OPTIONS] 

		--target  <FILE>  Target region file in bed format (required)

		--input   <FILE>  A list of bam files list. eg. bam_list_example.txt (required)

		--output  <FILE>  Directory for RPKM files (required)
		

2. Merge single sample RPKM files to a union data matrix

		usage: python erds_pe.py merge_rpkm [OPTIONS] 

		--rpkm_dir  <FILE>  Giving the RPKM files directory for taking data (required)

		--target   <FILE>  Target region file in bed format (required)

		--output  <FILE>  Output the data matrix file (optional)
		
		
3. Normalization

	This command is to normalize RPKM data matrix using principal component analysis.
	
		usage: python erds_pe.py svd [OPTIONS] 

		--rpkm_matrix  <FILE>  Giving the RPKM files PCA normalization (required)

		--output  <FILE>  Output the normalized data matrix file (optional)
		
2. Calling CNVs 

	This command is to call CNVs from pooled whole-exome sequencing samples.

		usage: python erds_pe.py discover [OPTIONS]  

		--params  <FILE>  Parameters files for HMM (required)

		--datafile  <FILE>  Normalized data matrix file (required)
		
		--output  <FILE>  Output the normalized data matrix file (optional)
		
		--sample  <STRING>  Giving a specific sample for calling (optional)

		--vcf  <FILE>  Taking SNV information from vcf file (optional)
		
			--hetsnp  <BOOL>  Using or not take heterogenous SNV information into HMM (optional, default FALSE)
		
			--tagsnp  <BOOL>  Using or not take tagSNP-copy number polymorphism information into HMM (optional, default FALSE)

		--tagsnp_file  <FILE>  A file records the linkage disequilibrium information between tagSNP and copy number polymorphism (optional)


#File Instruction

1. bam list file (three columns) 

	Column 1: path of .bam file

	Example: 

		/path/Sample1.bam

		/path/Sample2.bam

		/path/Sample3.bam
		...

2. pedigree file

	See (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml)

	Note: The Individual ID must be same as the @RG SM tag of the bam file.


#Contact
<yzhuangliu@gmail.com>
