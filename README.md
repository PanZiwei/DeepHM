## DeepH&M

README file created by Yu He on Oct 18, 2019

Please contact yu.he@wustl.edu with questions or for clarifications.


This folder "DeepHM" contains the code for running DeepH&M in the server. Download the code, which will create DeepHM folder .


a. After downloading the code in the folder, download mm9 datasets (for processing genomic features) and sample datasets (MeDIP-seq, MRE-seq, hmC-Seal, TAB-seq and WGBS) and MethPipe tool at the link http://wangftp.wustl.edu/~yhe/DeepHM_datasets/dataset.tar.gz. Uncompress the compressed .tar.gz file in DeepHM folder using command "tar -xzvf dataset.tar.gz". This will generate three folders: mm9, data, MethPipe. Folder mm9 lists files for processing genomic features for mm9. Folder data lists sample data for MeDIP-seq, MRE-seq, hmC-Seal, TAB-seq and WGBS data. Folder MethPipe lists the tool for running mlml. Below list detailed descriptions for files in folder mm9 and data.

1.folder mm9. The default species is mm9. If doing training for a different species, you need to provide following files for that species.

cpgIsland.bed: cpg island coordinates. Each column lists the coordinate for a cpg island (chr, start, end)

hmc_interval: hmc intervals used for dividing cpgs into multiple windows for data balance

mm9.fa: mm9 genome sequence

cpg_no_chrM: cpg coordinates (no ChrM). Each column lists the coordinate for a cpg (chr, start, end)

gc5Base.sort.bedGraph: gc percent. Each column lists the coordinate and gc perecent for the interval (chr, start, end, gc_percent)

mm9_chrom_sizes: chromosome size

ACGT_sites.bed: coordinates for ACGT restriction sites

GCGC_sites.bed: coordinates for GCGC restriction sites

CCGC_sites.bed: coordinates for CCGC restriction sites

CCGG_sites.bed: coordinates for CCGG restriction sites


2.folder data. All files are sorted.

hmC-Seal.bedGraph: bedGraph file for hmc-seal data. Each colum lists the coordinate and number of reads overlapping the interval

MeDIP.bedGraph: bedGraph file for medip-seq data

MRE.bed: bed file for mre-seq data. Each column is a read

TAB_data: tab-seq data. Each colum lists coordinates, 5hmC level, coverage

WGBS_data: wgbs data. Each colum lists coordinates, total methylation level, coverage




b. Dependencies for DeepHM include BedTools, GSL(v2.2.1), bash, python3 and R. Please make sure these are installed and working before trying to run DeepHM.




c. After finishing the setup, follow the steps below to run DeepH&M. 

Run following scripts in DeepHM folder by first setting PWD=`pwd`.




1: process genomic features for mm9 or others

Script: 

```bat
bash process_genome.sh mm9_folder genome_process_folder
```

Inputs:

mm9_folder: genome folder which stores mm9 data

genome_process_folder: output folder where processed mm9 data will be stored




2: process tab-seq, wgbs, medip,hmcSeal and mre data.

The process data step will use ``libgsl.so.19``in **GSL shared libraries** (https://www.gnu.org/software/gsl/doc/html/usage.html#shared-libraries), you need to locate the file ``libgsl.so.19``:

2.1 Installing Anaconda and GSL(v2.2.1) package

2.2 Searched for ``libgsl.so.19`` and found the location: ``/path/to/library``

2.3 Run the following scripts to adds the path of the library to an environment variable. This step is only valid for the current session. It will not work for other users. It will not work once you log off and on again.

``bat
LD_LIBRARY_PATH=$LD_LIBRARY_PATH: /path/to/library
export LD_LIBRARY_PATH
``

2.4 Run the process data script:

```bat
bash process_methylation_data.sh mm9_folder data_folder data_process_folder genome_process_folder
```

Inputs:

mm9_folder: genome folder which stores mm9 data

data_folder: data folder which stores tab-seq, wgbs, medip,hmcSeal and mre data

data_process_folder: output folder where processed data will be stored  

genome_process_folder: processed genome folder which stores output data from step 1


3: select high cov data, normalize and balance data .

Script: 
```bat
bash tune_data.sh data_process_folder data_tune_folder mm9_folder
```

Inputs:

data_process_folder: processed data folder which stores output data from step 2

data_tune_folder: output folder where tuned data will be stored  

mm9_folder: genome folder which stores mm9 data




4: train cpg and dna model .

Script: 

```bat
bash train_cpg_dna_model.sh data_tune_folder train_folder genome_process_folder
```

Inputs:

data_tune_folder: tuned data folder which stores output data from step 3

train_folder: output folder where training data will be stored  

genome_process_folder: processed genome folder which stores output data from step 1



5: train joint model .

Script: 

```bat
bash train_joint_model.sh train_folder cpg_model dna_model genome_process_folder
```

Inputs:

train_folder: training data folder which stores output data from step 4

cpg_model: trained cpg model from step 4

dna_model: trained dna model from step 4  

genome_process_folder: processed genome folder which stores output data from step 1



6: predict from cpg, dna and joint model . For predicting on current dataset, use following script.

Script: 

```bat
bash predict_model.sh train_folder cpg_model dna_model joint_model pred_folder data_process_folder data_tune_folder mm9_folder genome_process_folder
```

Inputs:

train_folder: training data folder which stores output data from step 4 

cpg_model: trained cpg model from step 4 

dna_model: trained cpg model from step 4 

joint_model: trained cpg model from step 5 

pred_folder: output folder where prediction data will be stored   

data_process_folder: processed data folder which stores output data from step 2 

data_tune_folder: tuned data folder which stores output data from step 3 

mm9_folder: genome folder which stores mm9 data

genome_process_folder: processed genome folder which stores output data from step 1

Outputs:
all_coord_pred_hmc_mc: predicted 5hmc and 5mc levels. Each column is chr, start, end, predicted 5hmc, predicted 5mc



7: predict from cpg, dna and joint model . For predicting on new dataset, use following scripts.

Scripts:

1).process data for new dataset. Input variables are same as above.

```bat
bash process_methylation_data.sh mm9_folder data_folder data_process_folder genome_process_folder
```

2).predict data for new dataset using trained model. Input variables are same as above.

```bat
bash predict_model.sh train_folder cpg_model dna_model joint_model pred_folder data_process_folder data_tune_folder mm9_folder genome_process_folder
```

d. An example to train DeepHM model using provided mm9 and data folder:

```bat
bash process_genome.sh mm9 genome_process_folder

bash process_methylation_data.sh mm9 data data_process_folder genome_process_folder

bash tune_data.sh data_process_folder data_tune_folder mm9

bash train_cpg_dna_model.sh data_tune_folder train_folder genome_process_folder

bash train_joint_model.sh train_folder cpg_model-epoch20 dna_model-epoch30 genome_process_folder

bash predict_model.sh train_folder cpg_model-epoch20 dna_model-epoch30 joint_model-epoch40 pred_folder data_process_folder data_tune_folder mm9 genome_process_folder
```
