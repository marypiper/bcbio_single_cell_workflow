# bcbio_singlecell_rnaseq workflow

1. Ask client for the following:
	- How many samples were sequenced?
	- What were the sample indices used?
	- How many cells were encapsulated and sequenced per sample?

2. Acquire data from sequencing core:

	- **Bauer sequencing core** uses Basespace and to download the sequencing files, follow the code below:

		```
		wget https://da1s119xsxmu0.cloudfront.net/sites/knowledgebase/API/08052014/Script/	BaseSpaceRunDownloader_v2.zip
		unzip BaseSpaceRunDownloader_v2.zip
		python BaseSpaceRunDownloader_v2.py -r <Run ID> -a <access token>
		```
		
		The option `-r` is the number in the basespace url and the `access token` is something you have to get for your basespace account. The files output will be BCL files that require demultiplexing with the `bcl2fastq` tool (instructions below).
		
		More information can be found at https://developer.basespace.illumina.com/docs/content/documentation/authentication/obtaining-access-tokens.

	- **DFCI sequencing center (Zach)** will output the FASTQ files (already demultiplexed).
	- **Biopolymers sequencing facility** will sometimes output BCL and sometimes FASTQ, so necessary to check the files
	- **Broad Institute** has their own single cell distribution platfrom

3. If downloaded sequencing files are BCL format, then need to convert to FASTQ. To do this log on to Orchestra or O2 to run `bcl2fastq`. Currently the module is only available on Orchestra.

	- Change directories to the sequencing folder downloaded from the facility. The folder should be arranged according to the image below for NextSeq or MiniSeq:
	

		<img src="../img/sequencing_dir_org.png" width="400">
	
		*Image acquired from [bcl2fastq documentation](../docs/bcl2fastq2_guide_15051736_v2.pdf).* 
	
	- Load `bcl2fastq` module and convert files to FASTQ by using the following command:
		
		```
		bcl2fastq \
		--use-bases-mask y*,y*,y*,y* \
		--mask-short-adapter-reads 0 \
		--minimum-trimmed-read-length 0
		```
		
		More information regarding the `bcl2fastq` command and directory structures for other sequencing machines can be found in the [documentation](../docs/bcl2fastq2_guide_15051736_v2.pdf). 
		
4. The output files should be in the `BaseCalls` directory. For each file of sequenced reads, there should be four associated FASTQ files (R1-R4). These 
	
	- R1 (61 bp Read 1): sequence of the read
	- R2 (8 bp Index Read 1 (i7)): cellular barcode - which cell read originated from
	- R3 (8 bp Index Read 2 (i5)): library index - which sample read originated from
	- R4 (14 bp Read 2): read 2 and barcode/UMI - remaining cellular barcode and UMI - which transcript read originated from (to find PCR duplicates)

	The reads for each sequence are depicted in the image below:

	<img src="../img/sc_seq_method.png" width="800">
	
5. To quickly count up the FASTQ entries in a file:

	```
	gzip -cd klein-v3_R3.fq.gz | head -40000 | awk 'NR % 4 == 2' | sort | uniq -c | awk 	'{ print $2 "," $1}' | sort -t"," -n --key=2 | tail -5
	```
	
	*NOTE: ``awk 'NR % 4 == 2'` gets every 4th line starting from the 2nd, which is a useful trick when you want to count up FASTQ file entries*
	
6. Use the information from the client to construct the metadata table according to the specifications detailed at [https://github.com/hbc/bcbioSingleCell](https://github.com/hbc/bcbioSingleCell).

6. 
