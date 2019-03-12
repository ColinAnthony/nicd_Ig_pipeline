# nicd_Ig_pipeline
This wrapper script was written to formalize/standardize the NICD antibody sequencing pipeline

This script depends on:
 
    linux
    python3.6 or later (https://www.anaconda.com/distribution/)
    python libraries:
        docopt   (conda install docopt)
        pathlib  (conda install pathlib)
        pandas   (conda install pandas)
    vsearch (https://github.com/torognes/vsearch)
    PEAR (https://cme.h-its.org/exelixis/web/software/pear/)
    SONAR scripts for annotating nAb lineages (available from xxx)

This script requires a strict naming format:

    STUDYPID_visit_wpi_chain_primername using
    
    "_" as a delimeter
    do not include "-" in your file name. 
    These will be converted to "_" and may cause the sctipts to break
    
      eg: CAP255_4180_80wpi_heavy_C5
      Where study = CAP
      PID = 255 (Zero padded to three digits - ie: 008 not 8)
      visit = 4180
      wpi = 080wpi (Zero padded to three digits - ie: 080wpi not 80wpi)
      chain = heavy
      primername = C5

Sample names must match the names in the settings file exactly, for  the first 5 fields (by "_" separation)
    
    ie: CAP255_4180_080wpi_heavy_C5_S2_L001_R1_001.fastq
        CAP255_4180_080wpi_heavy_C5_S2_L001_R2_001.fastq

The settings file file must contain these headings:
    (A template for the settings file is available in this repo)
    
    sample_name	
    sonar_1_version	
    lineage	
    primer_name	
    time_point	
    run_step1	
    run_step2	
    run_step3	
    known_mab_name

If you have two mAb sequences that you want to run a sample against on Sonar P2,
include an entry for each mAb, as shown in the settings template


**To run the pipeline:**
    
* create a project folder (usually named after the participant)
* inside this folder, create a folder called 0_new_data
    * copy your .zip archive containing the Illumina paired end files into this folder
        * also accepts .gz compressed files or .fastq files
    * prepare your settings.csv file, indicating which steps to run on which samples
        * if running additional samples, set previous samples to '0' in the three run_step columns
* create a folder called `mab_sequences` which contains a fasta file with all the mAb sequencse you will need for sonar P2
 * run the wrapper script:
    * it is recommened to use `screen` or `nohup` as the run times will be long.
    
    `sudo apt install screen`
    
    `screen -S <job_name>`
    
     `ig_pipeline.py -p <project_path> -s <settings_file>`
 
 * check the log file to see details of the processes that were run
 * check your output files
 * do your analysis
 
 