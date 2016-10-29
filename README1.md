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

		--params  <FILE>  Parameters file for HMM (required)

		--datafile  <FILE>  Normalized data matrix file (required)
		
		--output  <FILE>  Output the normalized data matrix file (optional)
		
		--sample  <STRING>  Giving a specific sample for calling (optional)

		--vcf  <FILE>  Taking SNV information from vcf file (optional)
		
			--hetsnp  <BOOL>  Using or not take heterogenous SNV information into HMM (optional, default FALSE)
		
			--tagsnp  <BOOL>  Using or not take tagSNP-copy number polymorphism information into HMM (optional, default FALSE)

		--tagsnp_file  <FILE>  A file records the linkage disequilibrium information between tagSNP and copy number polymorphism (optional)


#File Instruction

1. bam list file (three columns) 

	Column 1: Sample_name
	Column 2: Path of bam files
	Column 3: The population of the corresponding sample. e.g CEU, YRI, CHB etc.

	Example: 

		NA06984	/data/rjtan/1000GP/exome/NA06984.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam	CEU
		NA06985	/data/rjtan/1000GP/exome/NA06985.mapped.ILLUMINA.bwa.CEU.exome.20130415.bam	CEU
		NA06994	/data/rjtan/1000GP/exome/NA06994.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam	CEU
		NA07000	/data/rjtan/1000GP/exome/NA07000.mapped.ILLUMINA.bwa.CEU.exome.20130415.bam	CEU
		...

#Quick Start & a example

1. Getting RD and transformming to RPKM format

	This command is to extract the read depth (RD) signal from the BAM files and to calculate RPKM values for all samples. 

		usage: python erds_pe.py rpkm [OPTIONS] 

		--target  <FILE>  Target region file in bed format (required)

		--input   <FILE>  A list of bam files list. eg. bam_list_example.txt (required)

		--output  <FILE>  Directory for RPKM files (required)
		
1. Step 1:

	python erds_pe.py rpkm
	input $bamlist_file
	target $target_file
	output $rpkm_files
	
2. Step 2:

	python erds_pe.py merge_rpkm
	--rpkm_dir $rpkm_files
	--target $target_file
	
3. Step 3:

	python erds_pe.py svd
	--rpkm_matrix $RPKM_matrix.raw.filtered

4. Step 4:

	python erds_exome.py discover
	params params.txt
	--datafile $RPKM_matrix.raw.filtered.SVD
	--sample NA12878
	--vcf=$snv_vcf_file
	--hetsnp True
	--tagsnp Ture
	--tagsnp_file=$tagsnp_file
	--output NA12878.pooled.Het.Tag.cnv
  
#Contact
<renjie.tan@outlook.com>
