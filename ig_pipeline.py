#!/usr/local/bin/python3.6
"""This is a wrapper to automate the running of the NGS nAb processing and analysis pipeline on their cluster.

    Usage:
        ig_pipeline.py (-p <project_path>) (-s <settings_file>) [-f <fasta_file>] [-r]
        ig_pipeline.py -h
        ig_pipeline.py -v

    Options:
        -p --path               The path to the project folder, where the folders will be created
        -s --settings           The path and file name of the settings csv file
        -f --fasta_file         The path and name of the fasta file with your mAb sequences
        -r --run_sonar2_trunc   Will run the Sonar P2 trunc call
                                    (Use this flag only if you have double peaks in your divergence plots)
        -v --version            Show the script version number
        -h --help               Show this screen.
"""
import sys
import os
import subprocess
import collections
import datetime
import random
import string
from itertools import groupby
import pathlib
import filecmp
# external libraries
from docopt import docopt
import pandas as pd


__author__ = 'Colin Anthony'


class dircmp(filecmp.dircmp):
    """
    https://stackoverflow.com/questions/4187564/recursive-dircmp-compare-two-
    directories-to-ensure-they-have-the-same-files-and
    Compare the content of dir1 and dir2. In contrast with filecmp.dircmp, this
    subclass compares the content of files with the same path.
    """
    def phase3(self):
        """
        Find out differences between common files.
        Ensure we are using content comparison with shallow=False.
        """
        fcomp = filecmp.cmpfiles(self.left, self.right, self.common_files, shallow=False)
        self.same_files, self.diff_files, self.funny_files = fcomp


def is_same(dir1, dir2):
    """
    https://stackoverflow.com/questions/4187564/recursive-dircmp-compare-two-
    directories-to-ensure-they-have-the-same-files-and
    Compare two directory trees content.
    :param dir1: (str) path to dir 1
    :param dir2: (str) path to dir 2
    Return False if they differ, True is they are the same.
    """
    compared = dircmp(dir1, dir2)
    if (compared.left_only or compared.right_only or compared.diff_files
        or compared.funny_files):
        return False
    for subdir in compared.common_dirs:
        if not is_same(os.path.join(dir1, subdir), os.path.join(dir2, subdir)):
            return False
    return True


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for k, v in my_gen:
        v = v.replace("-", "")
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = v.upper()

    return dct


def disk_space_checker(path):
    # check for free disk space
    disk_space_object = os.statvfs(path)
    gb = (1024 * 1024) * 1024
    # always keep a bit of free space for the safety reasons
    safety_buffer = 10
    total_space = disk_space_object.f_blocks * disk_space_object.f_frsize / gb - safety_buffer
    free_space = disk_space_object.f_bfree * disk_space_object.f_frsize / gb - safety_buffer
    percent_free = round((free_space / total_space) * 100, 2)
    return free_space, percent_free


def folder_size_checker(start_path='.'):
    """
    function to get the size of a folder on disk, taken
    from https://stackoverflow.com/questions/1392413/calculating-a-directorys-size-using-python
    :param start_path:
    :return:
    """
    total_size = sum(os.path.getsize(f) for f in os.listdir(start_path) if os.path.isfile(f))

    return total_size


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
        job_list_entry = sample_id, [lineage, time_point, chain], [run_step1, run_step2, run_step3], primer_name, known_mab_name
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
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "crd3", "work").mkdir(mode=0o777, parents=True,
                                                                                         exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "crd3", "output").mkdir(mode=0o777, parents=True,
                                                                                           exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "fullmab", "work").mkdir(mode=0o777, parents=True,
                                                                                            exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, known_mab_name, "fullab", "output").mkdir(mode=0o777, parents=True,
                                                                                             exist_ok=True)


def unzip_files(path, logfile):
    """
    a function to unpack .zip archives and decompress any .gz files
    :param path: the path to the project folder
    :param logfile: (str) path and name of the logfile
    :return:
    """
    # Bites to Gb adjustment
    gb = (1024 * 1024) * 1024

    new_data = pathlib.Path(path, "0_new_data")
    if not new_data.is_dir():
        print(f"'0_new_data' folder not found\n"
              f"You need to create a folder '0_new_data' in the project folder, ie:\n{new_data}"
              f"\n and copy your data in there")
        print("making the folder for you")
        new_data.mkdir(mode=0o777, parents=True, exist_ok=True)
        sys.exit("exiting")

    # find zip files
    search_zip = list(new_data.glob("*.zip"))
    if search_zip:
        for file in search_zip:
            # get for free disk space
            free_space, percent_free = disk_space_checker(path)
            # get file size
            zip_size = os.path.getsize(file) / gb
            if free_space < zip_size * 5:
                print(f"not enough free space to unzip\nFree space is {free_space}Gb")
                with open(logfile, "a") as handle:
                    handle.write(f"# Not enough free space to unzip\nFree space is {free_space}Gb {percent_free}%\n")
                sys.exit("exiting")
            else:
                print("archived file detected")
                cmd_unzip = f"unzip {file} -d {new_data}"
                try:
                    subprocess.call(cmd_unzip, shell=True)
                    with open(logfile, "a") as handle:
                        handle.write("# unzipping archive\n")
                        handle.write(cmd_unzip + "\n")
                        os.unlink(str(file))
                        # log command
                        handle.write("# removing original archive\n")
                except subprocess.CalledProcessError as e:
                    print(e)
                    print("unzipping archive failed")
                    if list(pathlib.Path(path, "0_new_data").glob("*.gz")):
                        continue
                    else:
                        sys.exit("exiting")

    # find .gz files
    search_gz = list(new_data.glob("*.gz"))

    # get needed space to unpack gz files
    total_un_gz_size = 0
    # get free space
    free_space, percent_free = disk_space_checker(path)
    for file in search_gz:
        total_un_gz_size += round(int(subprocess.check_output(f"zcat {file} | wc --bytes", shell=True)) / gb, 3)
    if free_space < total_un_gz_size * 2:
        print("not enough free space to unzip\n")
        with open(logfile, "a") as handle:
            handle.write(f"# Not enough free space to all the gz files\nNeed {total_un_gz_size} just to unzip"
                         f"\nFree space is {free_space}")
        sys.exit("exiting")

    # unpack the gz files
    if search_gz:
        for file in search_gz:
            # do replace all '-' with '_' in file name
            cmd_rename = f"rename -v '-' '_' {file}"
            os.chmod(str(file), 0o777)
            subprocess.call(cmd_rename, shell=True)
            # log command
            with open(logfile, "a") as handle:
                handle.write("# replacing '-' with '_'\n")
                handle.write(cmd_rename + "\n")
        # set job id
        gunzip_job_name = 'gunzip'
        run_gunzip = pathlib.Path(path, "scripts", "run_gunzip.sh")

        with open(run_gunzip, "w") as handle:
            handle.write("#!/bin/sh\n")
            handle.write("#SBATCH -J gzip\n")
            handle.write("#SBATCH --mem=1000\n\n")
            handle.write(f"gunzip {str(new_data)}/*.gz")

        os.chmod(run_gunzip, 0o777)
        cmd_gunzip = f"sbatch -J {gunzip_job_name} {run_gunzip}"
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
        # primer_name = job_entry[3]
        known_mab_name = job_entry[4]
        run_steps = job_entry[2]
        chain = job_entry[1][2]
        dir_with_raw_files = pathlib.Path(path, job_entry[1][0], job_entry[1][1], job_entry[1][2], "1_raw_data")
        dir_with_sonar1_files = pathlib.Path(path, job_entry[1][0], job_entry[1][1], job_entry[1][2], "4_dereplicated")
        dir_with_sonar2_files = pathlib.Path(path, job_entry[1][0], job_entry[1][1], job_entry[1][2], "4_dereplicated",
                                             f"5_{known_mab_name}", "output")
        if 1 in run_steps:
            for i, step in enumerate(run_steps):
                if step == 1:
                    if i == 0:
                        command_call_processing.append([sample_name, dir_with_raw_files])
                        search_raw_path = list(dir_with_raw_files.glob("*.fastq*"))
                        if search_raw_path:
                            n += 1

                    elif i == 1:
                        command_call_sonar_1.append([sample_name, dir_with_sonar1_files, chain])
                        search_sonar1_path = list(dir_with_sonar1_files.glob("*.fasta*"))
                        if search_sonar1_path:
                            n += 1
                    elif i == 2:
                        command_call_sonar_2.append([sample_name, dir_with_sonar1_files, chain, known_mab_name])
                        search_sonar2_path = list(dir_with_sonar2_files.glob("*.fasta"))
                        if search_sonar2_path:
                            n += 1
                    else:
                        print("this should not happen\nnumber or run steps larger than expected\n")
                        sys.exit("exiting")
        else:
            continue

    # exit if no files were found
    if n == 0:
        print("No files found in target directories")
        check_bad_person = list(path.glob("*.fastq*"))
        if len(check_bad_person) > 0:
            print("raw data found in project folder but not in '0_new_data' folder\nmove files to '0_new_data'")
        else:
            sys.exit("exiting")
    else:
        with open(logfile, "a") as handle:
            handle.write(f"# files found in target folder\n")
        print("files found")

    return command_call_processing, command_call_sonar_1, command_call_sonar_2


def step_1_run_sample_processing(path, command_call_processing, logfile):
    """
    function to process the raw fastq files into dereplicated fasta files for sonar1
    :param path: (str) path to the project folder
    :param command_call_processing: (list of lists) each list contains the sample name and path to the raw data
    :param logfile: (str) path and name of the logfile
    :return: none
    """
    # Bites to Gb adjustment
    gb = (1024 * 1024) * 1024

    for item in command_call_processing:
        sample_name = item[0]
        with open(logfile, "a") as handle:
            handle.write(f"running Sonar P1 on {sample_name}\n")

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
                # get for free disk space
                free_space, percent_free = disk_space_checker(path)
                zip_size = os.path.getsize(file) / gb
                if free_space < zip_size * 2:
                    print(f"not enough free space to gunzip\nFree space is {free_space}Gb")
                    with open(logfile, "a") as handle:
                        handle.write(
                            f"# Not enough free space to gunzip\nFree space is {free_space}Gb {percent_free}%\n")
                    sys.exit("exiting")
                else:
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
        if not search_raw_files:
            print(f"No fastq files in target folder for {sample_name}\nTrying next sample\n")
            continue

        for file_r1 in search_raw_files:
            file_r2 = str(file_r1).replace("_R1.fastq", "_R2.fastq")
            file_r2 = pathlib.Path(file_r2)
            # if can't find the matched paired R2 file, skip and go to next R1 file
            if not file_r2.is_file():
                print(f"R2 paired file not found\nSkipping this sample: {file_r1}")
                continue
            else:
                # get for free disk space
                free_space, percent_free = disk_space_checker(path)
                file_size = os.path.getsize(file_r1) / gb * 2
                if free_space < file_size * 4:
                    print(f"not enough disk space to process the raw files\nNeed{file_size * 4}Gb\n"
                          f"you have{free_space}Gb free")
                    with open(logfile, "a") as handle:
                        handle.write(
                            f"# not enough disk space to process the raw files\nNeed {file_size * 4}Gb\n"
                            f"you have{free_space}Gb free\n")
                    sys.exit("exiting")

                # if both R1 and R2 files are present, continue with the processing
                print("running PEAR")
                pear_outfile = file_r1.name.replace("_R1.fastq", "")
                pear_outfile = f"{pear_outfile}"
                merged_file_name = pathlib.Path(merged_folder, pear_outfile)

                # create SLURM file to run PEAR on the cluster
                id_prefix = sample_name[3:6]
                unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
                # set job id
                pear_job_name = f"{id_prefix}{unique_suffix}PR"
                run_pear = pathlib.Path(dir_with_raw_files, "run_pear.sh")
                with open(run_pear, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=node01\n")
                    handle.write("#SBATCH --mem=1000\n\n")
                    handle.write(f"/opt/conda2/pkgs/pear-0.9.6-2/bin/pear -f {file_r1} -r {file_r2} "
                                 f"-o {merged_file_name} -p 0.001 -j 4 -q 20 -n 300")
                os.chmod(str(run_pear), 0o777)
                with open(logfile, "a") as handle:
                    handle.write(f"# running PEAR command from file:\n{run_pear}\n")
                cmd_pear = f"sbatch -J {pear_job_name} {run_pear}"
                subprocess.call(cmd_pear, shell=True)
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
                    cmd_zip = f"gzip {file_r1} {file_r2}"
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
                run_fastq_fasta = pathlib.Path(merged_folder, "run_convert_to_fasta.sh")
                with open(run_fastq_fasta, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=node01\n")
                    handle.write("#SBATCH --mem=1000\n\n")
                    handle.write(f"/opt/conda2/pkgs/vsearch-2.4.3-0/bin/vsearch --fastq_filter {fastq} "
                                 f"--fastaout  {fasta} --fasta_width 0 --notrunclabels --threads 4")
                os.chmod(run_fastq_fasta, 0o777)
                with open(logfile, "a") as handle:
                    handle.write(f"# running fastq_to_fasta command from file:\n{run_fastq_fasta}\n")
                cmd_fastq_fasta = f"sbatch -J {fastq_fasta_job_name} {run_fastq_fasta}"
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

    # move fastas to target dir and do dereplication if needed
    with open(logfile, "a") as handle:
        handle.write(f"# dereplicating files\n")

    # collect all the files that will be dereplicated
    files_to_derep_dict = collections.defaultdict(list)
    for item in command_call_processing:
        sample_name = item[0]
        dir_with_raw_files = item[1]
        parent_path = dir_with_raw_files.parent
        files_to_derep_dict[str(parent_path)].append(sample_name)

    # for each target folder, asign the relevant variables
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

        # check if you need to concatenate files
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
        # set dereplicated outfile name
        dereplicated_file = pathlib.Path(derep_folder, f"{fasta_to_derep_name_stem}_unique.fasta")
        # create fastq to fasta SLURM file
        run_derep = pathlib.Path(fasta_folder, "run_derep.sh")
        with open(run_derep, "w") as handle:
            handle.write("#!/bin/sh\n")
            handle.write("##SBATCH -w, --nodelist=node01\n")
            handle.write("#SBATCH --mem=1000\n\n")
            handle.write(f"/opt/conda2/pkgs/vsearch-2.4.3-0/bin/vsearch --sizeout --derep_fulllength "
                         f"{file_to_dereplicate} --output {dereplicated_file} --fasta_width 0 "
                         f"--notrunclabels --threads 4 ")
        os.chmod(run_derep, 0o777)
        with open(logfile, "a") as handle:
            handle.write(f"# running dereplication command from file:\n{str(file_to_dereplicate)}\n")
        cmd_derep = f"sbatch -J {derep_job_name} {run_derep}"

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


def step_2_run_sonar_p1(command_call_sonar_1, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_1: list of samples to run sonar1 on
    :param logfile: (str) path and name of the logfile
    :return:
    """
    gb = (1024 * 1024) * 1024

    for item in command_call_sonar_1:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        project_folder = dir_with_sonar1_files.parents[-2]
        print(project_folder)
        input("enter")
        dir_with_sonar1_work = pathlib.Path(dir_with_sonar1_files, "work")
        dir_with_sonar1_output = pathlib.Path(dir_with_sonar1_files, "output")

        with open(logfile, "a") as handle:
            handle.write(f"running Sonar P1 on {sample_name}\n")

        # check whether there is already sonar P1 output in the target folders (if this is a rerun)
        search_work = list(dir_with_sonar1_work.glob("*"))
        search_output = list(dir_with_sonar1_output.glob("*"))
        if search_work:
            with open(logfile, "a") as handle:
                handle.write(f"# removing existing files from sonar1 target directory {dir_with_sonar1_work}\n")
            for file in search_work:
                os.unlink(str(file))
        if search_output:
            with open(logfile, "a") as handle:
                handle.write(f"# removing existing files from sonar1 target directory {dir_with_sonar1_output}\n")
            for file in search_output:
                os.unlink(str(file))

        sonar_version = item[2]
        project_folder = dir_with_sonar1_files.parents[3]

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
        else:
            # check that you have enough space
            file = search_derep_fastas[0]
            # get for free disk space
            free_space, percent_free = disk_space_checker(project_folder)
            fasta_size = os.path.getsize(file) / gb
            if free_space < fasta_size * 10:
                print(f"not enough free space to run sonar P1\nFree space is {free_space}Gb"
                      f"# you expected to need up to {fasta_size * 10}Gb for sonar P1")
                with open(logfile, "a") as handle:
                    handle.write(
                        f"# Not enough free space to sonar P1\nFree space is {free_space}Gb {percent_free}%\n"
                        f"# you expected to need up to {fasta_size * 10}Gb for sonar P1")
                sys.exit("exiting")

        id_prefix = sample_name[3:6]
        unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
        # set job id
        sonar1_job_name = f"{id_prefix}{unique_suffix}S1"

        # get sonar P1 version specific commands
        sonar_settings = {"heavy": ["H", "-lib /opt/conda2/pkgs/sonar/germDB/IgHJ.fa "
                                         "-dlib /opt/conda2/pkgs/sonar/germDB/IgHD.fa "
                                         "-clib /opt/conda2/pkgs/sonar/germDB/IgHC_CH1.fa -callFinal"],
                          "kappa": ["K", "-lib /opt/conda2/pkgs/sonar/germDB/IgKJ.fa -noD -noC -callFinal"],
                          "lambda": ["L", "-lib /opt/conda2/pkgs/sonar/germDB/IgLJ.fa -noD -noC -callFinal"]}
        sonar_version_setting = sonar_settings[sonar_version]

        run_sonar_p1 = pathlib.Path(project_folder, "scripts", f"run_sonar_P1_{sonar_version}.sh")
        with open(run_sonar_p1, "w") as handle:
            handle.write("#!/bin/sh\n")
            handle.write("#SBATCH -w, --nodelist=bio-linux\n")
            handle.write("#SBATCH --mem=4000\n\n")
            handle.write(f"python2 /opt/conda2/pkgs/sonar/annotate/1.1-blast_V.py -locus {sonar_version_setting[0]} "
                         f"-callJ -jArgs {sonar_version_setting[1]}\n")
            handle.write(f"sbatch -J {sonar1_job_name} {run_sonar_p1}")
        os.chmod(str(run_sonar_p1), 0o777)
        with open(logfile, "a") as handle:
            handle.write(f"# running Sonar1 command from file:\n{str(run_sonar_p1)}\n")
        try:
            # change into sonar P1 targer directory
            os.chdir(dir_with_sonar1_files)
            # subprocess.call(run_sonar, shell=True)
        except subprocess.CalledProcessError as e:
            print(e)
            print("Sonar P1 encountered an error\ntrying next sample")
            with open(logfile, "a") as handle:
                handle.write(f"Sonar P1 encountered \n{e}\n")
            continue
        # return cwd to the project folder
        os.chdir(project_folder)


def step_3_run_sonar_2(command_call_sonar_2, fasta_sequences, run_sonar2_trunc, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_2: list of samples to run sonar2 on
    ":param fasta_sequences: (dict) dict of all the mAb fasta sequences (full and crd3)
    :param run_sonar2_trunc: (Bool) True if flag is set, default is False, will run sonar_trunc call if True
    :param fasta_sequences: (dict) dict containing key = seq name, value = sequence for all known mAb sequences
    :param logfile: (str) path and name of the logfile
    :return:
    """
    gb = (1024 * 1024) * 1024

    for item in command_call_sonar_2:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        parent_dir = dir_with_sonar1_files.parents[3]
        chain = item[2]
        known_mab_name = item[3]
        target_folder_full_ab = pathlib.Path(f"5_{known_mab_name}", "full_ab")
        target_folder_crdh3 = pathlib.Path(f"5_{known_mab_name}", "crd3")
        fullab_name = known_mab_name + f"_{chain}_fullab"
        cdr3_name = known_mab_name + f"_{chain}_cdr3"
        mab_sequence_fullab = fasta_sequences[fullab_name]
        mab_sequence_cdr3 = fasta_sequences[cdr3_name]

        # check for mab_sequences folder
        fasta_sequences = pathlib.Path(parent_dir, "mab_sequences")
        if not fasta_sequences.is_dir():
            print(f"'mab_sequences' folder not found\n"
                  f"You need to create a folder 'mab_sequences' in the project folder, ie:\n{fasta_sequences}"
                  f"\nand copy your fasta file in there")
            print("making the folder for you")
            fasta_sequences.mkdir(mode=0o777, parents=True, exist_ok=True)
            sys.exit("exiting")

        # check for Sonar P1 files
        dir_with_sonar1_work = pathlib.Path(dir_with_sonar1_files, "work")
        dir_with_sonar1_output = pathlib.Path(dir_with_sonar1_files, "output")
        for sonar1_dir in [dir_with_sonar1_work, dir_with_sonar1_output]:
            sonar_p1_search = sonar1_dir.glob("*")
            if not sonar_p1_search:
                print(f"No Sonar P1 files detected in: {dir_with_sonar1_files}")
                sys.exit()

        with open(logfile, "a") as handle:
            handle.write(f"running Sonar P1 on {sample_name}\n")
        to_run = [(target_folder_full_ab, fullab_name, mab_sequence_fullab),
                  (target_folder_crdh3, cdr3_name, mab_sequence_cdr3)]
        mab = ["fullab", "cdr3"]
        for i, (target_folder, mab_name, mab_seq) in enumerate(to_run):
            # make fasta file for fullab run
            mab_name_file = f"{item}.fasta"
            mab_name_file = pathlib.Path(target_folder, mab_name_file)
            with open(mab_name_file, 'w') as handle:
                handle.write(f">{mab_name}\n{mab_sequence_fullab}\n")

            # make Sonar P2 script
            if run_sonar2_trunc:
                sonar_p2_run = f"{known_mab_name}_{mab[i]}_sonar_p2_run_trunc.sh"
                sonar_p2_run = pathlib.Path(parent_dir, sonar_p2_run)
                with open(sonar_p2_run, 'w') as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=bio-linux\n")
                    handle.write("#SBATCH --mem=4000\n\n")
                    handle.write(f"perl /opt/conda2/pkgs/sonar/lineage/2.1-calculate_id-div.pl "
                                 f"-a {mab_name_file} -g /opt/conda2/pkgs/sonar/germDB/IgHKLV_cysTruncated.fa "
                                 f"-ap muscle -t 4\n")
            else:
                sonar_p2_run = f"{known_mab_name}_{mab[i]}_sonar_p2_run.sh"
                with open(sonar_p2_run, 'w') as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=bio-linux\n")
                    handle.write("#SBATCH --mem=4000\n\n")
                    handle.write(f"perl /opt/conda2/pkgs/sonar/lineage/2.1-calculate_id-div.pl "
                                 f"-a {mab_name_file} -ap muscle -t 4\n")

            # check that target dirs are empty
            dir_with_sonar2_work = pathlib.Path(target_folder, "work")
            dir_with_sonar2_output = pathlib.Path(target_folder, "output")
            # check that target dirs have unmodified sonar P1 output
            sonar_p1_alread_cp_work = is_same(dir_with_sonar1_work, dir_with_sonar2_work)
            sonar_p1_alread_cp_output = is_same(dir_with_sonar1_output, dir_with_sonar2_output)
            copy_sonar_p1 = True
            for sonar2_dir in [dir_with_sonar2_work, dir_with_sonar2_output]:
                sonar_p2_search = sonar2_dir.glob("*")
                if sonar_p2_search:
                    if sonar_p1_alread_cp_work and sonar_p1_alread_cp_output:
                        print("unmodified sonar P1 output found in sonar P2 targer folder\nNot re-copying the files")
                        copy_sonar_p1 = False
                    else:
                        print(f"Modified sonar P1 output detected in sonar P2 target folder {item}\ndeleting files")
                        for file in sonar_p2_search:
                            os.unlink(str(file))

                else:
                    # if sonar P1 output not already (unmodified) in target dir
                    if copy_sonar_p1:
                        # check for free disk space
                        free_space, percent_free = disk_space_checker(parent_dir)
                        sonar_p1_size = round(folder_size_checker(dir_with_sonar1_output) / gb, 3)
                        if free_space < sonar_p1_size * 3:
                            print(f"not enough free space to run sonar P2\nFree space is {free_space}Gb"
                                  f"# you expected to need up to {sonar_p1_size * 3}Gb for sonar P2")
                            with open(logfile, "a") as handle:
                                handle.write(f"# Not enough free space to sonar P2\n"
                                             f"Free space is {free_space}Gb {percent_free}%\n"
                                             f"# you expected to need up to {sonar_p1_size * 3}Gb for sonar P2")
                            sys.exit("exiting")
                        # copy files to desired directory
                        cmd_cope_work = f"cp {dir_with_sonar1_work} {target_folder}"
                        cmd_cope_output = f"cp {dir_with_sonar1_output} {target_folder}"
                        try:
                            subprocess.call(cmd_cope_work, shell=True)
                            subprocess.call(cmd_cope_output, shell=True)
                        except subprocess.CalledProcessError as e:
                            print(e)
                            print("There was an error copying the sonar P1 files to the sonar P2 directory")
                            with open(logfile, "a") as handle:
                                handle.write(f"Error copying Sonar P1 output to Sonar P2 folder\n{e}\n")
                            continue

                        sonar_p1_alread_cp_work = is_same(dir_with_sonar1_work, dir_with_sonar2_work)
                        sonar_p1_alread_cp_output = is_same(dir_with_sonar1_output, dir_with_sonar2_output)
                        if not sonar_p1_alread_cp_output and sonar_p1_alread_cp_work:
                            print("Copied files are different from sonar P1 original files\nCopy must have failed")
                            sys.exit("exiting")
            # run sonar2
            os.chdir(target_folder)
            os.chmod(str(sonar_p2_run), 0o777)
            try:
                subprocess.call(sonar_p2_run, shell=True)
            except subprocess.CalledProcessError as e:
                print(e)
                print("There was an error running sonar P2\nTrying next sample")
                with open(logfile, "a") as handle:
                    handle.write(f"There was an error running sonar P2\n{e}\n")
                continue


def main(path, settings, fasta_file=None, run_sonar2_trunc=False):
    """
    create folder structure for a sequencing project
    :param path: (str) path to where the folders should be made
    :param settings: (str) the path and csv file with the run settings
    :param fasta_file: (str) the path and name of the fasta file with your mAb sequences
    :param run_sonar2_trunc: (Bool) set flag to true to run sonar_trunc
    """
    path = pathlib.Path(path).absolute()
    # change the cwd to project path
    os.chdir(path)

    # check disk space
    free_space, percent_free = disk_space_checker(path)
    print(f"Free disk space remaining: {free_space}Gb ({percent_free}%)")

    # set paths and initialize the log file
    project_name = path.parts[-1]
    settings = pathlib.Path(settings).absolute()
    settings_dataframe = pd.read_csv(settings, sep=None, engine='python')
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(path, f"{time_stamp}_log_file.txt")
    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for project: {project_name}\n")
        handle.write(f"Free space at start = {free_space}Gb ({percent_free}%)\n\n")
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
        step_1_run_sample_processing(path, command_call_processing, log_file)

    # only run sonar1 if one or more files were specified
    if command_call_sonar_1:
        # get rid of duplicate entries in sonar1 list, if present
        dedup_sonar1_call = []
        removed_sonor1_duplicates_check = []
        for item in command_call_sonar_1:
            sample_name = "_".join(item[0].split("_")[:-1])
            sonar1_dir = item[1]
            chain = item[2]
            if sonar1_dir in removed_sonor1_duplicates_check:
                continue
            else:
                removed_sonor1_duplicates_check.append(sonar1_dir)
                dedup_sonar1_call.append([sample_name, sonar1_dir, chain])

        step_2_run_sonar_p1(dedup_sonar1_call, log_file)

    # only run sonar 2 if one or more files were specified
    if command_call_sonar_2:
        if fasta_file is not None:
            fasta_sequences = fasta_to_dct(fasta_file)
            step_3_run_sonar_2(command_call_sonar_2, fasta_sequences, run_sonar2_trunc, log_file)
        else:
            print("no fasta file specified for Sonar P2 to run\n")
            sys.exit("exiting")
        pass

    print("Done")


if __name__ == "__main__":
    args = docopt(__doc__, version='V0.1.0')

    main(path=args["<project_path>"], settings=args["<settings_file>"], fasta_file=args["<fasta_file>"],
         run_sonar2_trunc=args["--run_sonar2_trunc"])
