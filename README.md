# PipelineProject_Cian_Dotson
This repository is a walkthrough of an analysis of the transcriptomes of individuals 2 and 6 days post infection (dpi) with Huamn herpevirus, or Huamn cytomegalovirus, or HCMV.

# **Step 1: Downloading the Raw Fastq Files**
The raw fastq files for RNA-Seq data we are interested in are found in these four links:

Donor 1 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

Donor 1 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

Donor 2 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

Donor 2 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

To get the SRA Normalized data for each of these accessions, simply use `wget`:

`wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030`

`wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033`

`wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044`  

`wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045`

These files are not yet readable as fastq files, but you can extract the fastq files by using `fasterq-dump` NCBI's SRA-toolkit:

`fasterq-dump SRR5660030`

`fasterq-dump SRR5660033`

`fasterq-dump SRR5660044`

`fasterq-dump SRR5660045`
For ease, a subset of the data are stored in the ./reference directory when you download this direcetory.

# **Step 2: Downloading the Dependencies**:
In order to run this pipeline, you need the following software installed. The links to their installation instructions are listed in order of appearance in the pipeline below:

**python version 3.10.12**: https://www.python.org/downloads/ 

Downloading this version of python will have the modules **argparse**, **sys**, **os**, **re**, **csv**, and **statitsics** already installed

**biopython version 1.83**: https://biopython.org/wiki/Download 

**kallisto version 0.51.1**: https://github.com/pachterlab/kallisto 

**R version 4.4.2**: https://repo.miserver.it.umich.edu/cran/ 

**sleuth version 0.30.1**: https://github.com/pachterlab/sleuth

**dplyr version 1.1.4**: https://cran.r-project.org/web/packages/dplyr/readme/README.html

**bowtie2 version 2.4.4**: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2

**SPAdes version 4.0.0**: https://github.com/ablab/spades/releases/tag/v4.1.0

**blast+ version 2.12.0**: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# **Step 3: Using demo_transcriptomics.py**:
To run 




  
