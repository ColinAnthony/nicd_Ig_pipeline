#!/anaconda/bin/python3
"""This is a wrapper to automate the NGS nAb processing and analysis pipeline.

    Usage:
         ig_pipeline.py (-p <project_path>) (-s <settings_file>)
         ig_pipeline.py -h
         ig_pipeline.py -v

    Options:
         -p --path          The path to the project folder, where the folders will be created
         -s --settings      The path and file name of the settings csv file
         -v --version       Show the script version number
         -h --help          Show this screen.
"""
# builtin libraries
import sys
import os
import subprocess
import collections
import datetime
import random
import string
# external libraries
from docopt import docopt
import pathlib
import pandas as pd


__author__ = 'Colin Anthony'


def settings_checker(settings_dataframe):
    """
    function to check the format of the settings csv file
    :param settings_dataframe: dataframe object of settings csv file
    :return: dict of settings dataframe, key = index for each sample entry
    """
    # check that the settings file had the correct headings
    headings = list(settings_dataframe)
    expected_headings = ["sample_name", "sonar_1_version", "lineage", "primer_name", "time_point", "run_step1",
                         "run_step2", "run_step3", "known_mab_name"]
    for item in expected_headings:
        if item not in headings:
            print(f"{item} heading not found in settings file.\n\nPlease fix the settings file before running"
                  f"\nExpected headings are: {expected_headings}")
            sys.exit("exiting")
    # stop nan type being set as float
    for header in ["sample_name", "sonar_1_version", "lineage", "primer_name", "time_point", "known_mab_name"]:
            settings_dataframe[header] = settings_dataframe[header].astype(str)

    # convert dataframe to dict
    settings_dict = settings_dataframe.to_dict(orient='index', into=dict)
    run_steps = []
    for job_entry, job_settings in settings_dict.items():
        sample_id = job_settings["sample_name"]
        parts = sample_id.split("_")
        if len(parts) != 5:
            print(f"your sample name was {sample_id}\n")
            print("Sample name was not correctly formatted\n"
                  "Use '_' as the field delimeter\n"
                  "Sample name must be: pid_visit_wpi_chain_primername\n"
                  "eg: CAP255_4180_80wpi_heavy_C5")
            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[0]) != 6:
            print(f"your sample name was {sample_id}\n The PID must be zero padded: eg: CAP008, not CAP8")
            print("Sample name was not correctly formatted\n"
                  "Use '_' as the field delimeter\n"
                  "Sample name must be: pid_visit_wpi_chain_primername\n"
                  "eg: CAP255_4180_80wpi_heavy_C5")
            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[1]) != 4:
            print(f"your sample name was {sample_id}\n The visit code must be in the 2000 format, not P2V0")
            print("Sample name was not correctly formatted\n"
                  "Use '_' as the field delimeter\n"
                  "Sample name must be: pid_visit_wpi_chain_primername\n"
                  "eg: CAP255_4180_80wpi_heavy_C5")
            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[2]) != 6 and parts[2][-3:].lower() != "wpi":
            print(f"your sample name was {sample_id}\n The wpi code must be zero padded and end in wpi eg: 080wpi")
            print("Sample name was not correctly formatted\n"
                  "Use '_' as the field delimeter\n"
                  "Sample name must be: pid_visit_wpi_chain_primername\n"
                  "eg: CAP255_4180_80wpi_heavy_C5")
            sys.exit("exiting")

        chain = job_settings["sonar_1_version"].lower()
        chain_options = ["heavy", "kappa", "lambda"]
        if chain not in chain_options:
            print(f"sonar_1_version {chain} was not one of the expected options: {chain_options}")
            sys.exit("exiting")
        lineage = job_settings["lineage"].lower()
        if lineage == "nan":
            print(f" lineage must be specified")
            sys.exit("exiting")
        time_point = job_settings["time_point"].lower()
        if time_point == "nan":
            print(f" time_point must be specified")
            sys.exit("exiting")
        primer_name = job_settings["primer_name"].upper()
        if primer_name == "NAN":
            print(f" primer_name must be specified")
            sys.exit("exiting")
        known_mab_name = job_settings["known_mab_name"].upper()
        if known_mab_name == "NAN":
            print(f" known_mab_name must be specified")
            sys.exit("exiting")
        run_steps.append(job_settings["run_step1"])
        run_steps.append(job_settings["run_step2"])
        run_steps.append(job_settings["run_step3"])

    if 1 not in set(run_steps):
        print("no run options were set\nSet either run_step1, run_step2, or run_step3 to 1, for at least one entry")

    return settings_dict


def extract_settings_make_folders(path, settings_dict):
    """
    function to extract the settings information for each entry and call the make folders function
    :param path: path to project folder
    :param settings_dict: the settings in dictionary format
    :return: list of lists for the job entries
    """
    list_all_jobs_to_run = []
    for job_entry, job_settings in settings_dict.items():
        sample_id = job_settings["sample_name"]
        known_mab_name = job_settings["known_mab_name"].upper()
        run_step1 = job_settings["run_step1"]
        run_step2 = job_settings["run_step2"]
        run_step3 = job_settings["run_step3"]
        lineage = job_settings["lineage"].lower()
        time_point = job_settings["time_point"].lower()
        chain = job_settings["sonar_1_version"].lower()
        primer_name = job_settings["primer_name"].upper()
        job_list_entry = sample_id, [lineage, time_point, chain], [run_step1, run_step2, run_step3], \
                         primer_name, known_mab_name
        list_all_jobs_to_run.append(job_list_entry)

        # make the folders if they don't already exist
        step_0_make_folders(path, lineage, time_point, chain, known_mab_name)

    return list_all_jobs_to_run


def step_0_make_folders(path, lineage, time_point, chain, known_mab_name):
    """
    function to create the nested folder structure for the nAb pileline, if they don't already exist
    :param path: (str) the path to the project folder
    :param lineage: (str) Ab lineage name
    :param time_point: (str) the time point the NGS data was sampled from
    :param chain: (str) the Ab chain (heavy, kappa, lambda)
    :param known_mab_name: (str) the primer name used for target amplification
    :return: None
    """

    known_mab_name = "5_" + known_mab_name
    pathlib.Path(path, "scripts").mkdir(mode=0o777, parents=True, exist_ok=True)

    pathlib.Path(path, lineage, time_point, chain, "1_raw_data").mkdir(mode=0o777, parents=True, exist_ok=True)

    pathlib.Path(path, lineage, time_point, chain, "2_merged_filtered").mkdir(mode=0o777, parents=True, exist_ok=True)

    pathlib.Path(path, lineage, time_point, chain, "3_fasta").mkdir(mode=0o777, parents=True, exist_ok=True)

    pathlib.Path(path, lineage, time_point, chain, "4_dereplicated", "work").mkdir(mode=0o777, parents=True,
                                                                                   exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "4_dereplicated", "output").mkdir(mode=0o777, parents=True,
                                                                                     exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "crdh3", "work").mkdir(mode=0o777, parents=True,
                                                                                          exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "crdh3", "output").mkdir(mode=0o777, parents=True,
                                                                                            exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "full_ab", "work").mkdir(mode=0o777, parents=True,
                                                                                            exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "full_ab", "output").mkdir(mode=0o777, parents=True,
                                                                                              exist_ok=True)


def unzip_files(path, logfile):
    """
    a function to unpack .zip archives and decompress any .gz files
    :param path: the path to the project folder
    :param logfile: (str) path and name of the logfile
    :return:
    """
    new_data = pathlib.Path(path, "0_new_data")
    search_zip = list(new_data.glob("*.zip"))
    if search_zip:
        for file in search_zip:
            if file.suffix == ".zip":
                print("archived file detected")
                cmd_jar = f"unzip {file} -d {new_data}"
                print(cmd_jar)
                try:
                    subprocess.call(cmd_jar, shell=True)
                    with open(logfile, "a") as handle:
                        handle.write("# unzipping archive\n")
                        handle.write(cmd_jar + "\n")
                        # os.unlink(str(file))
                        # log command
                        handle.write("# removing original archive\n")
                except subprocess.CalledProcessError as e:
                    print(e)
                    print("unzipping archive failed")
                    if list(pathlib.Path(path, "0_new_data").glob("*.gz")):
                        continue
                    else:
                        sys.exit("exiting")

    search_gz = list(new_data.glob("*.gz"))
    if search_gz:
        for file in search_gz:
            # do replace all '-' with '_' in file name
            # Todo: change to nicd version of rename
            # cmd = f"rename -v '-' '_' {file}"
            cmd_rename = f"rename 's/-/_/g' {file}"
            os.chmod(str(file), 0o777)
            subprocess.call(cmd_rename, shell=True)
            # log command
            with open(logfile, "a") as handle:
                handle.write("# replacing '-' with '_'\n")
                handle.write(cmd_rename + "\n")

        # Todo: make job name
        gunzip_job_name = ''
        run_gunzip = pathlib.Path(path, "scripts", "run_gunzip.sh")

        with open(run_gunzip, "w") as handle:
            # handle.write("#!/bin/sh\n")
            # handle.write("#SBATCH -J gzip\n")
            # handle.write("#SBATCH --mem=1000\n\n")
            handle.write(f"gunzip {str(new_data)}*.gz")
        # Todo: change to nicd version
        os.chmod(run_gunzip, 0o777)
        # cmd_gunzip = f"sbatch -J {gunzip_job_name} {run_gunzip}"
        cmd_gunzip = f"{run_gunzip}"
        subprocess.call(cmd_gunzip, shell=True)
        # log command
        with open(logfile, "a") as handle:
            handle.write("# uncompressing .gz files\n")
            handle.write(cmd_gunzip + "\n")


def move_raw_data(path, settings_dataframe, logfile):
    """
    function to move raw data to the correct folder within the project
    :param path: path to the project folder
    :param settings_dataframe: the dataframe object form the settings file
    :param logfile: (str) path and name of the logfile
    :return: a list of all the files that were moved
    """
    raw_files = pathlib.Path(path, "0_new_data").glob("*.fastq")

    # rename raw files
    for n, file in enumerate(raw_files):
        new_name = str(file).replace("-", "_")
        os.rename(str(file), new_name)

    # find where to move each raw file and move it to the right folder
    raw_files = list(pathlib.Path(path, "0_new_data").glob("*.fastq"))
    if not raw_files:
        print("No fastq (fastq.gz/fastq.zip) files in 0new_data")
        with open(logfile, "a") as handle:
            handle.write("# No fastq (fastq.gz/fastq.zip) files in 0new_data\n")

    for file in raw_files:
        full_name = file
        name = file.stem
        parts = name.split("_")
        search_name = "_".join(parts[:5])
        copy_name = "_".join(parts[:5]) + f"_{parts[7]}.fastq"
        fields = settings_dataframe.loc[settings_dataframe['sample_name'] == search_name].head(1)
        if fields.empty:
            print(f"Your sample name {search_name}\nwas not found in the settings 'sample_name' column\n"
                  f"Please fix the file name or settings file accordingly")
            with open(logfile, "a") as handle:
                handle.write(f"Your sample name {search_name}\nwas not found in the settings 'sample_name' column\n")
                handle.write(f"Please fix the file name or settings file accordingly\nexiting")
            sys.exit("exiting")

        lineage = fields.iloc[0]["lineage"].lower()
        chain = fields.iloc[0]["sonar_1_version"].lower()
        time_point = fields.iloc[0]["time_point"].lower()
        destination = pathlib.Path(path, lineage, time_point, chain, "1_raw_data", copy_name)

        # copy the file to the correct folder
        if pathlib.Path(destination).is_file():
            print(f"The file \n{full_name}\nis already in destination folder")
            check = True
            while check:
                answer = input(f"Type: 'yes' to overwrite the existing file in the destination folder'\n"
                               f"Type: 'no' to remove the file from 0new_data: ").lower()
                if answer not in ["yes", "no"]:
                    print("\nResponse was not a valid answer, enter 'yes', or 'no'\n")
                    pass
                elif answer == "yes":
                    cmd1 = f"mv {full_name} {destination}"
                    subprocess.call(cmd1, shell=True)
                    check = False
                else:
                    check = False
                    os.unlink(str(full_name))

        else:
            print("moving file")
            print(full_name, "\n")

            cmd = f"mv {full_name} {destination}"
            subprocess.call(cmd, shell=True)
            with open(logfile, "a") as handle:
                handle.write(f"# moving files to target folder\n{cmd}\n")


def make_job_lists(path, list_all_jobs_to_run, logfile):
    """
    function to generate a list of samples to run for each step in pipeline
    :param path: (pathlib object) the project folder path
    :param list_all_jobs_to_run:
    :param logfile: (str) path and name of the logfile
    :return: (list) of jobs to run for each step
    """
    # counter for checking if there are files to run the pipeline on
    n = 0
    # get job list entries
    command_call_processing = []
    command_call_sonar_1 = []
    command_call_sonar_2 = []

    for job_entry in list_all_jobs_to_run:
        sample_name = job_entry[0]
        dir_with_raw_files = pathlib.Path(path, job_entry[1][0], job_entry[1][1], job_entry[1][2], "1_raw_data")
        dir_with_sonar1_files = pathlib.Path(path, job_entry[1][0], job_entry[1][1], job_entry[1][2], "4_dereplicated")

        chain = job_entry[1][2]
        primer_name = job_entry[3]
        known_mab_name = job_entry[4]
        run_steps = job_entry[2]
        if 1 in run_steps:
            for i, step in enumerate(run_steps):
                if step == 1:
                    if i == 0:
                        command_call_processing.append([sample_name, dir_with_raw_files])
                    elif i == 1:
                        command_call_sonar_1.append([sample_name, dir_with_sonar1_files, chain])
                    elif i == 2:
                        command_call_sonar_2.append([sample_name, dir_with_sonar1_files, chain, primer_name,
                                                     known_mab_name])
                    else:
                        print("this should not happen\nnumber or run steps larger than expected\n")
                        sys.exit("exiting")
        else:
            continue

        # check that files exist in target folder
        search_raw_path = list(dir_with_raw_files.glob("*.fastq*"))
        search_sonar1_path = list(dir_with_sonar1_files.glob("*.fasta*"))

        # increment n with the count of files found
        n += len(search_raw_path)
        n += len(search_sonar1_path)

    # exit of no files were found
    print("Checking whether files are in target folder")
    if n == 0:
        print("No files found in target directories")
        sys.exit("exiting")
    else:
        with open(logfile, "a") as handle:
            handle.write(f"# files found in target folder\n")
        print("files found")

    return command_call_processing, command_call_sonar_1, command_call_sonar_2


def step_1_run_sample_processing(command_call_processing, logfile):
    """
    function to process the raw fastq files into dereplicated fasta files for sonar1
    :param command_call_processing: (list of lists) each list contains the sample name and path to the raw data
    :param logfile: (str) path and name of the logfile
    :return: none
    """
    for item in command_call_processing:
        sample_name = item[0]
        with open(logfile, "a") as handle:
            handle.write(f"working on {sample_name}\n")

        print(f"processing {sample_name}")

        dir_with_raw_files = item[1]
        parent_path = dir_with_raw_files.parent
        merged_folder = pathlib.Path(parent_path, "2_merged_filtered")
        fasta_folder = pathlib.Path(parent_path, "3_fasta")
        derep_folder = pathlib.Path(parent_path, "4_dereplicated")

        # check for zipped files
        search_zip_raw_files = list(dir_with_raw_files.glob(f"{sample_name}_R*.fastq.gz"))
        if search_zip_raw_files:
            for file in search_zip_raw_files:
                cmd_gunzip = f"gunzip {file}"
                subprocess.call(cmd_gunzip, shell=True)

            # remove old merged file
            print(f"removing old merged file for {sample_name}")
            for file in merged_folder.glob(f"{sample_name}*.fastq*"):
                os.unlink(str(file))

            # remove old fasta file
            print(f"removing old fasta file for {sample_name}")
            for file in fasta_folder.glob(f"{sample_name}*.fasta*"):
                os.unlink(str(file))

            # remove old dereplicated file
            print(f"removing old dereplicated file for {sample_name}")
            trunc_name = "_".join(sample_name.split("_")[:-1])
            for file in derep_folder.glob(f"{trunc_name}*_unique.fasta"):
                os.unlink(str(file))

        search_raw_files = list(dir_with_raw_files.glob(f"{sample_name}_R1.fastq"))
        for file_R1 in search_raw_files:
            file_R2 = str(file_R1).replace("_R1.fastq", "_R2.fastq")
            file_R2 = pathlib.Path(file_R2)
            # if can't find the matched paired R2 file, skip and go to next R1 file
            if not file_R2.is_file():
                print(f"R2 paired file not found\nSkipping this sample: {file_R1}")
                continue
            else:
                # if both R1 and R2 files are present, continue with the processing
                print("running PEAR")
                pear_outfile = file_R1.name.replace("_R1.fastq", "")
                pear_outfile = f"{pear_outfile}"
                merged_file_name = pathlib.Path(merged_folder, pear_outfile)

                # create SLURM file to run PEAR on the cluster
                # Todo: change to nicd version
                id_prefix = sample_name[3:6]
                unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
                # set job id
                pear_job_name = f"{id_prefix}{unique_suffix}PR"
                run_pear = pathlib.Path(dir_with_raw_files, "run_pear.sh")
                with open(run_pear, "w") as handle:
                    # handle.write("#!/bin/sh\n")
                    # handle.write("##SBATCH -w, --nodelist=node01\n")
                    # handle.write("#SBATCH --mem=1000\n\n")
                    handle.write(f"pear -f {file_R1} -r {file_R2} -o {merged_file_name} -p 0.001 -j 4 -q 20 -n 300 -u 5")
                os.chmod(str(run_pear), 0o777)
                with open(logfile, "a") as handle:
                    handle.write(f"# running PEAR command from file:\n{run_pear}\n")
                # cmd_pear = f"sbatch -J {pear_job_name} {run_pear}"
                cmd_pear = f"{run_pear}"
                # subprocess.call(cmd_pear, shell=True)
                try:
                    pear_output = subprocess.check_output(cmd_pear, shell=True)
                    pear_output = pear_output.decode("utf-8")
                    with open(logfile, "a") as handle:
                        handle.write(pear_output)
                        handle.write("\n")
                except subprocess.CalledProcessError as e:
                    print(e)
                    print("PEAR encountered an error\ntrying next sample")
                    with open(logfile, "a") as handle:
                        handle.write(f"PEAR failed\n{e}\n")
                    continue

                # compress 1_raw_data if the merging was successful
                merged_files_list = list(pathlib.Path(merged_folder).glob(f"{sample_name}.assembled.fastq"))
                if merged_files_list:
                    cmd_zip = f"gzip {file_R1} {file_R2}"
                    subprocess.call(cmd_zip, shell=True)

                # remove unmerged files
                unassmebled_search = pathlib.Path(merged_folder).glob(f"{sample_name}*unassembled*.fastq")
                discarded_search = pathlib.Path(merged_folder).glob(f"{sample_name}*discarded.fastq")
                for file in unassmebled_search:
                    os.unlink(str(file))

                for file in discarded_search:
                    os.unlink(str(file))

                # convert to fasta
                if len(merged_files_list) > 1:
                    # this should not be possible, since any existing file should have been overwritten
                    print(f"multiple merged files found for this sample:\n{merged_files_list}\n")
                    sys.exit("exiting")
                elif len(merged_files_list) == 0:
                    print(f"no merged files found for this sample:\n{sample_name}")
                    sys.exit("exiting")
                else:
                    fastq = merged_files_list[0]
                fasta = f"{str(fastq.stem)}.fasta"
                fasta = pathlib.Path(fasta_folder, fasta)

                # create fastq to fasta SLURM file
                id_prefix = sample_name[3:6]
                unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
                # set job id
                fastq_fasta_job_name = f"{id_prefix}{unique_suffix}Fa"
                # Todo: change to nicd version
                run_fastq_fasta = pathlib.Path(merged_folder, "run_convert_to_fasta.sh")
                with open(run_fastq_fasta, "w") as handle:
                    # handle.write("#!/bin/sh\n")
                    # handle.write("##SBATCH -w, --nodelist=node01\n")
                    # handle.write("#SBATCH --mem=1000\n\n")
                    handle.write(f"vsearch --fastq_filter {fastq} --fastaout  {fasta} --fasta_width 0 --notrunclabels "
                                 f"--threads 4")
                os.chmod(run_fastq_fasta, 0o777)
                with open(logfile, "a") as handle:
                    handle.write(f"# running fastq_to_fasta command from file:\n{run_fastq_fasta}\n")
                # cmd_fastq_fasta = f"sbatch -J {fastq_fasta_job_name} {run_fastq_fasta}"
                cmd_fastq_fasta = f"{run_fastq_fasta}"
                try:
                    subprocess.call(cmd_fastq_fasta, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                except subprocess.CalledProcessError as e:
                    print(e)
                    print("vsearch fastq to fasta encountered an error\ntrying next sample")
                    with open(logfile, "a") as handle:
                        handle.write(f"vsearch fastq to fasta failed\n{e}\n")
                    continue
                # compress pear if the merged file was successfully converted to a fasta file
                if fasta.is_file():
                    # add file to list of fastas to derep (or just be moved to derep folder if only a single fasta)
                    cmd_zip_merged = f"gzip {fastq}"
                    subprocess.call(cmd_zip_merged, shell=True)
                else:
                    print("could not find the fasta file\nconversion of fastq to fasta might have failed\n")
                    sys.exit("exiting")

    # move fastas to target dir and do dereplication in needed
    with open(logfile, "a") as handle:
        handle.write(f"# dereplicating files\n")

    files_to_derep_dict = collections.defaultdict(list)
    for item in command_call_processing:
        sample_name = item[0]
        dir_with_raw_files = item[1]
        parent_path = dir_with_raw_files.parent
        files_to_derep_dict[str(parent_path)].append(sample_name)

    for folder, sample_name_list in files_to_derep_dict.items():
        fasta_folder = pathlib.Path(folder, "3_fasta")
        derep_folder = pathlib.Path(folder, "4_dereplicated")

        search_fasta_folder = list(pathlib.Path(fasta_folder).glob("*.fasta"))
        primers = []
        name_stem = ''
        for i, fasta_file in enumerate(search_fasta_folder):
            name = fasta_file.stem.replace(".assembled", "").split("_")
            if i == 0:
                name_stem += "_".join(name[:-1])
            primer = name[-1]
            primers.append(primer)
        primers_code = "_".join(primers)

        search_fasta_folder = list(pathlib.Path(fasta_folder).glob("*.fasta"))
        fasta_to_derep_name_stem = f"{name_stem}_{primers_code}"
        if len(search_fasta_folder) > 1:
            concated_outfile = pathlib.Path(fasta_folder, f"{fasta_to_derep_name_stem}_concatenated.fasta")
            if concated_outfile.is_file():
                # remove existing file so that you don't accidentally concatenate in duplicate
                print(f"{concated_outfile}\nalready exists\nthis file will be overwritten")
                os.unlink(concated_outfile)

            # concatenate multiple fasta files
            for file in search_fasta_folder:
                concat_cmd = f"cat {file} >> {concated_outfile}"
                subprocess.call(concat_cmd, shell=True)
            file_to_dereplicate = concated_outfile
        # only one fasta file for this sample
        elif len(search_fasta_folder) == 1:
            # copy fasta to 4_dereplicated
            file_to_dereplicate = search_fasta_folder[0]
        else:
            print(f"no Fasta files were found\nskipping this sample: {sample_name_list}")
            continue

        # dereplicate sequences using vsearch
        id_prefix = name_stem[3:6]
        unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
        # set job id
        derep_job_name = f"{id_prefix}{unique_suffix}De"
        # Todo: change to nicd version
        # set dereplicated outfile name
        dereplicated_file = pathlib.Path(derep_folder, f"{fasta_to_derep_name_stem}_unique.fasta")
        # create fastq to fasta SLURM file
        run_derep = pathlib.Path(fasta_folder, "run_derep.sh")
        with open(run_derep, "w") as handle:
            # handle.write("#!/bin/sh\n")
            # handle.write("##SBATCH -w, --nodelist=node01\n")
            # handle.write("#SBATCH --mem=1000\n\n")
            handle.write(f"vsearch --sizeout --derep_fulllength {file_to_dereplicate} --output {dereplicated_file} "
                         f"--fasta_width 0 --notrunclabels --threads 4 ")
        os.chmod(run_derep, 0o777)
        with open(logfile, "a") as handle:
            handle.write(f"# running dereplication command from file:\n{str(file_to_dereplicate)}\n")
        # cmd_derep = f"sbatch -J {derep_job_name} {run_derep}"
        cmd_derep = f"{run_derep}"
        try:
            derep_output = subprocess.check_output(cmd_derep, shell=True)
            derep_output = derep_output.decode("utf-8")
            with open(logfile, "a") as handle:
                handle.write(f"{derep_output}\n")
        except subprocess.CalledProcessError as e:
            print(e)
            print("vsearch dereplication encountered an error\ntrying next sample")
            with open(logfile, "a") as handle:
                handle.write(f"vsearch dereplication failed\n{e}\n")
            continue
        # remove non-dereplicated file if you had to concatenate multiple files
        for file in pathlib.Path(fasta_folder).glob(f"*_concatenated.fasta"):
            os.unlink(str(file))


def step_2_run_sonar_1(command_call_sonar_1, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_1: list of samples to run sonar1 on
    :param logfile: (str) path and name of the logfile
    :return:
    """
    for item in command_call_sonar_1:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        parent = pathlib.Path(dir_with_sonar1_files).parent
        print(parent)
        input("enter")
        sonar_version = item[2]
        with open(logfile, "a") as handle:
            handle.write(f"# running sonar1 on {sample_name}\n")

        # check that only one file in target dir
        search_derep_fastas = list(dir_with_sonar1_files.glob("*.fasta"))
        if len(search_derep_fastas) > 1:
            print(f"multiple fasta files found in {dir_with_sonar1_files}\n"
                  f"Sonar1 can only run if there is one fasta file in this folder\ntrying next sample")
            with open(logfile, "a") as handle:
                handle.write(f"# multiple fasta files in  {dir_with_sonar1_files}\n")
            continue
        id_prefix = sample_name[3:6]
        unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
        # set job id
        sonar1_job_name = f"{id_prefix}{unique_suffix}S1"
        # get sonar P1 version specific commands
        sonar_settings = {"heavy": "some settings", "kappa": "some settings", "lambda": "some settings"}
        sonar_version_setting = sonar_settings[sonar_version]

        run_sonar = pathlib.Path(dir_with_sonar1_files, f"run_sonar_P1_{sonar_version}.sh")
        with open(run_sonar, "w") as handle:
            # handle.write("#!/bin/sh\n")
            # handle.write("##SBATCH -w, --nodelist=node01\n")
            # handle.write("#SBATCH --mem=1000\n\n")
            # handle.write(sonar_version_setting)
            handle.write(f"sbatch -J {sonar1_job_name} {run_sonar}")
        os.chmod(str(run_sonar), 0o777)
        with open(logfile, "a") as handle:
            handle.write(f"# running Sonar1 command from file:\n{str(run_sonar)}\n")
        subprocess.call(run_sonar, shell=True)


def step_3_run_sonar_2(command_call_sonar_2, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_2: list of samples to run sonar2 on
    :param logfile: (str) path and name of the logfile
    :return:
    """
    for item in command_call_sonar_2:
        sample_name = item[0]
        dir_with_sonar2_files = item[1]
        sonar_version = item[2]
        primer_name = item[3]
        known_mab_name = item[4]
        known_mab_name_file = ''
        "#!/bin/sh"
        "##SBATCH -w, --nodelist=bio-linux"
        "#SBATCH --mem=26000"
        "perl /opt/conda2/pkgs/sonar/lineage/2.1-calculate_id-div.pl"
        f"-a {known_mab_name_file} -ap muscle"

        # check that target dirs are empty

        # check that there is enough disk space, send warning if not

        # copy files to desired directory

        # unzip files

        # run sonar2


def main(path, settings):
    """
    create folder structure for a sequencing project
    :param path: (str) path to where the folders should be made
    :param settings: (str) the path and csv file with the run settings
    """
    path = pathlib.Path(path).absolute()
    project_name = path.parent
    settings = pathlib.Path(settings).absolute()
    settings_dataframe = pd.read_csv(settings, sep=None, engine='python')
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(path, f"{time_stamp}_log_file.txt")
    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for {project_name}\n")
    os.chmod(str(log_file), 0o777)

    # check that settings file has the correct headings
    settings_dict = settings_checker(settings_dataframe)

    # make a list of the jobs to run and the settings required and make the required folders
    list_all_settings = extract_settings_make_folders(path, settings_dict)

    if not list_all_settings:
        print("No jobs were found in the settings file\nCheck that the file format was correct "
              "and that it contains entries")
        sys.exit("exiting")
    else:
        with open(log_file, "a") as handle:
            handle.write("# No fastq (fastq.gz/fastq.zip) files in 0new_data\n")

    # generate the job lists
    command_call_processing, command_call_sonar_1, command_call_sonar_2 = make_job_lists(path, list_all_settings,
                                                                                          log_file)
    # only run sample processing if one or more files were specified
    if command_call_processing:
        print("processing samples")

        # unzip_files(dir_with_files)
        unzip_files(path, log_file)

        # move files from 0_new_data, into correct directory, if necessary
        move_raw_data(path, settings_dataframe, log_file)
        # get rid of duplicate entries in sample processing list, if present
        command_call_processing = list(set(tuple(x) for x in command_call_processing))
        step_1_run_sample_processing(command_call_processing, log_file)

    # only run sonar1 if one or more files were specified
    if command_call_sonar_1:
        # get rid of duplicate entries in sonar1 list, if present
        dedup_sonar1_call = []
        removed_sonor1_duplicates_check = []
        for item in command_call_sonar_1:
            sample_name = "_".join(item[0].split("_")[:-1])
            sonar1_dir = item[1]
            chain = item[0]
            if sonar1_dir in removed_sonor1_duplicates_check:
                continue
            else:
                removed_sonor1_duplicates_check.append(sonar1_dir)
                dedup_sonar1_call.append([sample_name, sonar1_dir, chain])
    #
    #     step_2_run_sonar_1(command_call_sonar_1)
    #
    # # only run sonar 2 if one or more files were specified
    # if command_call_sonar_2:
    #     step_3_run_sonar_2(command_call_sonar_2)

    print("Done")


if __name__ == "__main__":
    args = docopt(__doc__, version='V0.1.0')

    main(path=args['<project_path>'], settings=args['<settings_file>'])
