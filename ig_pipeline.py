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
import struct
import time
import uuid
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
    if compared.left_only or compared.right_only or compared.diff_files or compared.funny_files:
        return False
    for subdir in compared.common_dirs:
        if not is_same(os.path.join(dir1, subdir), os.path.join(dir2, subdir)):
            return False
    return True


def read32(fileobj):
    """
    https://github.com/enthought/Python-2.7.3/blob/master/Lib/gzip.py
    :param fileobj: file object
    :return:
    """
    return struct.unpack("<I", fileobj.read(4))[0]


def get_uncompressed_size(file):
    """
    https://gist.github.com/ozanturksever/4968827
    :param file:
    :return:
    """
    with open(file, 'rb') as fileobj:
        fileobj.seek(-8, 2)
        isize = read32(fileobj)  # may exceed 2GB
    return isize


def disk_space_checker(path):
    """
    function to find how much disk space is free
    :param path: (str) path to check
    :return: (float) free space remaining in gb and percent of total
    """
    # check for free disk space
    disk_space_object = os.statvfs(path)
    gb = (1024 * 1024) * 1024
    # always keep a bit of free space for the safety reasons
    safety_buffer = 10
    total_space = disk_space_object.f_blocks * disk_space_object.f_frsize / gb - safety_buffer
    free_space = round(disk_space_object.f_bfree * disk_space_object.f_frsize / gb - safety_buffer, 3)
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


def settings_checker(settings_dataframe, logfile):
    """
    function to check the format of the settings csv file
    :param settings_dataframe: dataframe object of settings csv file
    :param logfile: (str) path and name of the logfile
    :return: dict of settings dataframe, key = index for each sample entry
    """
    # check that the settings file had the correct headings
    with open(logfile, "a") as handle:
        handle.write(f"\n# checking settings\n")
    headings = list(settings_dataframe)
    expected_headings = ["sample_name", "sonar_1_version", "lineage", "primer_name", "time_point", "run_step1",
                         "run_step2", "run_step3", "known_mab_name"]
    for item in expected_headings:
        if item not in headings:
            with open(logfile, "a") as handle:
                handle.write(f"# heading {item} not found in settings file.\n"
                             f"# Please fix the settings file before running\n"
                             f"# Expected headings are:\n\t{expected_headings}\n")
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
            with open(logfile, "a") as handle:
                handle.write(f"# Sample name was not correctly formatted: {sample_id}\n"
                             f"# expected eg: CAP255_4180_80wpi_heavy_C5\nUse '_' as the field delimeter\n")

            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[0]) != 6:
            with open(logfile, "a") as handle:
                handle.write(f"# Sample name was not correctly formatted: {sample_id}\n"
                             f"# expected the name to start with CAP\n"
                             f"# The PID must be zero padded: eg: CAP008, not CAP8\n"
                             f"# eg: CAP255_4180_80wpi_heavy_C5\n")
            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[1]) != 4:
            with open(logfile, "a") as handle:
                handle.write(f"# Sample name was not correctly formatted: {sample_id}\n"
                             f"# your sample name was {sample_id}\n"
                             f"#The visit code must be in the 2000 format, not P2V0\n"
                             f"# expected eg: CAP255_4180_80wpi_heavy_C5\n")
            sys.exit("exiting")
        elif parts[0][:3].upper() == "CAP" and len(parts[2]) != 6 and parts[2][-3:].lower() != "wpi":
            with open(logfile, "a") as handle:
                handle.write(f"# Sample name was not correctly formatted: {sample_id}\n"
                             f"# The wpi code must be zero padded and end in wpi eg: 080wpi"
                             f"# expected eg: CAP255_4180_80wpi_heavy_C5")
            sys.exit("exiting")

        chain = job_settings["sonar_1_version"].lower()
        chain_options = ["heavy", "kappa", "lambda"]
        if chain not in chain_options:
            with open(logfile, "a") as handle:
                handle.write(f"# Chain name was not in accepted list: {chain}\n"
                             f"expected: 'heavy', 'lambda' or 'kappa'")
            sys.exit("exiting")
        lineage = job_settings["lineage"].lower()
        if lineage == "nan":
            with open(logfile, "a") as handle:
                handle.write(f"# lineage variable not specified\n")
            sys.exit("exiting")
        time_point = job_settings["time_point"].lower()
        if time_point == "nan":
            with open(logfile, "a") as handle:
                handle.write(f"# time point variable not specified\n")
            sys.exit("exiting")
        primer_name = job_settings["primer_name"].upper()
        if primer_name == "NAN":
            with open(logfile, "a") as handle:
                handle.write(f"# primer_name variable not specified\n")
            sys.exit("exiting")
        known_mab_name = job_settings["known_mab_name"].upper()
        if known_mab_name == "NAN":
            with open(logfile, "a") as handle:
                handle.write(f"# known_mab_name variable not specified\n")
            sys.exit("exiting")
        run_steps.append(job_settings["run_step1"])
        run_steps.append(job_settings["run_step2"])
        run_steps.append(job_settings["run_step3"])

    if 1 not in set(run_steps):
        with open(logfile, "a") as handle:
            handle.write(f"# No run settings were detected\n"
                         f"# You must set run_step1, run_step2, or run_step3 to 1, for at least one entry")
        sys.exit("exiting")

    return settings_dict


def extract_settings_make_folders(path, settings_dict, logfile):
    """
    function to extract the settings information for each entry and call the make folders function
    :param path: path to project folder
    :param settings_dict: the settings in dictionary format
    :param logfile: (str) path and name of the logfile
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
        job_list_entry = sample_id, [lineage, time_point, chain], [run_step1, run_step2, run_step3], primer_name, \
                         known_mab_name
        list_all_jobs_to_run.append(job_list_entry)

        # make the folders if they don't already exist
        with open(logfile, "a") as handle:
            handle.write("\n# Making necessary folders for project\n")
        step_0_make_folders(path, lineage, time_point, chain, known_mab_name, logfile)

    return list_all_jobs_to_run


def step_0_make_folders(path, lineage, time_point, chain, known_mab_name, logfile):
    """
    function to create the nested folder structure for the nAb pileline, if they don't already exist
    :param path: (str) the path to the project folder
    :param lineage: (str) Ab lineage name
    :param time_point: (str) the time point the NGS data was sampled from
    :param chain: (str) the Ab chain (heavy, kappa, lambda)
    :param known_mab_name: (str) the primer name used for target amplification
    :param logfile: (str) path and name of the logfile
    :return: None
    """
    known_mab_name = "5_" + known_mab_name

    pathlib.Path(path, lineage, time_point, chain, "scripts").mkdir(mode=0o777, parents=True, exist_ok=True)

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

    new_data = pathlib.Path(path, "0_new_data")
    if not new_data.is_dir():
        with open(logfile, "a") as handle:
            handle.write("# 0_new_data' folder not found in the project folder\n#ie: {new_data}\n"
                         "# making the folder for you")
        new_data.mkdir(mode=0o777, parents=True, exist_ok=True)
        sys.exit("exiting")


def move_raw_data(path, settings_dataframe, logfile):
    """
    function to move raw data to the correct folder within the project
    :param path: path to the project folder
    :param settings_dataframe: the dataframe object form the settings file
    :param logfile: (str) path and name of the logfile
    :return: a list of all the files that were moved
    """

    raw_files = list(pathlib.Path(path, "0_new_data").glob("*.fastq*"))
    # rename raw files
    for file in raw_files:
        # fix file permissions for output
        os.chmod(str(file), 0o666)

        new_name = str(file).replace("-", "_")
        os.rename(str(file), new_name)

    # find where to move each raw file and move it to the right folder
    if not raw_files:
        with open(logfile, "a") as handle:
            handle.write("# No fastq/fastq.gz files in 0new_data\n")
    else:
        for file in raw_files:
            full_name = file
            name = file.stem
            parts = name.split("_")
            suff = file.suffixes
            suffix = ''.join(suff)
            if len(parts) < 6:
                with open(logfile, "a") as handle:
                    handle.write("# fastq filename not formatted correctly\n")
                sys.exit("exiting")

            if parts[-1] == 'R1' or parts[-1] == 'R2':
                # name already formatted
                copy_name = f"{name}{suffix}"
                search_name = "_".join(parts[:5])
            else:
                search_name = "_".join(parts[:5])
                copy_name = "_".join(parts[:5]) + f"_{parts[7]}{suffix}"

            fields = settings_dataframe.loc[settings_dataframe['sample_name'] == search_name].head(1)
            if fields.empty:
                with open(logfile, "a") as handle:
                    handle.write(f"# Your sample {search_name}\nwas not found in the settings 'sample_name' column\n")
                    handle.write(f"# Please fix the file name or settings file accordingly\n")
                sys.exit("exiting")

            lineage = fields.iloc[0]["lineage"].lower()
            chain = fields.iloc[0]["sonar_1_version"].lower()
            time_point = fields.iloc[0]["time_point"].lower()
            destination = pathlib.Path(path, lineage, time_point, chain, "1_raw_data", copy_name)

            # copy the file to the correct folder
            mv = f"mv {file} {destination}"
            with open(logfile, "a") as handle:
                handle.write(f"\n# moving files to target folder\n{mv}\n")

            try:
                with open(logfile, "a") as handle:
                    handle.write(f"# moving {full_name} to 1_raw_data failed\n{e}\n")
                subprocess.call(mv, shell=True)
            except subprocess.CalledProcessError as e:
                with open(logfile, "a") as handle:
                    handle.write(f"# moving {file} to 1_raw_data failed\n{e}\n")


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
                        with open(logfile, "a") as handle:
                            handle.write("# this should not happen\nnumber of run steps larger than expected\n")
                        sys.exit("exiting")
        else:
            continue

    # exit if no files were found
    if n == 0:
        check_bad_person = list(path.glob("*.fastq*"))
        check_good_person1 = list(pathlib.Path(path, "0_new_data").glob("*.fastq*"))
        if check_bad_person:
            with open(logfile, "a") as handle:
                handle.write("# raw data found in project folder but not in '0_new_data' folder\n"
                             "# move files to '0_new_data'")
            sys.exit("exiting")
        elif check_good_person1:
            pass
        else:
            with open(logfile, "a") as handle:
                handle.write("# No target files found")
            sys.exit("exiting")

    else:
        pass

    return command_call_processing, command_call_sonar_1, command_call_sonar_2


def raw_files_gz(chain_path, sample_name, dir_with_raw_files, script_folder, logfile):
    """

    :param chain_path:
    :param sample_name:
    :param dir_with_raw_files:
    :param script_folder:
    :param logfile:
    :return:
    """
    search_raw_files_fastq = list(dir_with_raw_files.glob(f"{sample_name}_R1.fastq"))
    if not search_raw_files_fastq:
        with open(logfile, "a") as handle:
            handle.write(f"# No fastq files in target folder for {sample_name}\n")
        raise IOError
    else:
        with open(logfile, "a") as handle:
            handle.write(f"# unzipped raw fastq file found in 1_raw_data for {sample_name}\n# zipping files")
        # set job id
        gz_unique_id = uuid.uuid4()
        gzip_job_name = f'gzipRaw'
        run_gzip = pathlib.Path(script_folder, "run_gzip_raw.sh")
        cmd_gz = f"gzip {dir_with_raw_files}/*.fastq"
        slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
        with open(run_gzip, "w") as handle:
            handle.write("#!/bin/sh\n")
            handle.write("#SBATCH -J gzip\n")
            handle.write("#SBATCH --mem=1000\n\n")
            handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
            handle.write(f"{cmd_gz}\n")
            handle.write(f"echo {gz_unique_id}")
        os.chmod(str(run_gzip), 0o777)

        with open(logfile, "a") as handle:
            handle.write(f"# gzip on raw file:\n{cmd_gz}\n")

        cmd_gzip = f"sbatch -J {gzip_job_name} {run_gzip} --parsable --wait"
        try:
            gzip_slurm_id = subprocess.check_output(cmd_gzip, shell=True).decode(sys.stdout.encoding).strip()
            gzip_slurm_id = gzip_slurm_id.split(" ")[-1]
            gz_slurm_out_file = pathlib.Path(chain_path, f"slurm-{gzip_slurm_id}.out")

        except subprocess.CalledProcessError as e:
            with open(logfile, "a") as handle:
                handle.write(f"# gzip on raw file failed\n{e}\n")
            raise
        search_raw_files = list(dir_with_raw_files.glob(f"{sample_name}_R1.fastq.gz"))

        if not search_raw_files:
            with open(logfile, "a") as handle:
                handle.write(f"# No fastq files in target folder for {sample_name}\n# Trying next sample\n")
            raise IOError

    return gz_slurm_out_file, gz_unique_id


def run_pear(chain_path, sample_name, script_folder, file_r1, file_r2, merged_file_name, itern, logfile):
    """

    :param chain_path:
    :param sample_name:
    :param script_folder:
    :param file_r1:
    :param file_r2:
    :param merged_file_name:
    :param itern:
    :param logfile:
    :return:
    """
    # create SLURM file to run PEAR on the cluster
    id_prefix = sample_name[3:6]
    unique_suffix = ''.join(random.choices(string.ascii_lowercase + string.digits, k=2)).lower()
    # set job id
    pear_unique_id = uuid.uuid4()
    pear_job_name = f"{id_prefix}{unique_suffix}PR"
    pear_script = pathlib.Path(script_folder, f"run_pear{str(itern)}.sh")
    pear = f"/opt/conda2/pkgs/pear-0.9.6-2/bin/pear  -f {file_r1} -r {file_r2} -o {merged_file_name} " \
        f"-p 0.001 -n 300\n"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(pear_script, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("##SBATCH -w, --nodelist=node01\n")
        handle.write("#SBATCH --mem=4000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{pear}\n")
        handle.write(f"echo {pear_unique_id}")
    os.chmod(str(pear_script), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# running PEAR\n{pear}\n\n")
    cmd_pear = f"sbatch -J {pear_job_name} {pear_script} --ntasks=1 --parsable --wait"
    try:
        pear_slurm_id = subprocess.check_output(cmd_pear, shell=True).decode(sys.stdout.encoding).strip()
        pear_slurm_id = pear_slurm_id.split(" ")[-1]
        pear_slurm_out_file = pathlib.Path(chain_path, f"slurm-{pear_slurm_id}.out")
    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"# pear failed\n{e}\n")
        raise

    return pear_slurm_out_file, pear_unique_id


def fastq2fasta(chain_path, sample_name, script_folder, fasta, merged_outfile, itern, logfile):
    """

    :param chain_path:
    :param sample_name:
    :param script_folder:
    :param fasta:
    :param merged_outfile:
    :param itern:
    :param logfile:
    :return:
    """
    # create fastq to fasta SLURM file
    id_prefix = sample_name[3:6]
    unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
    # set job id
    convert_unique_id = uuid.uuid4()
    fastq_fasta_job_name = f"{id_prefix}{unique_suffix}Fa"
    run_fastq_fasta = pathlib.Path(script_folder, f"run_convert_to_fasta{str(itern)}.sh")
    convert_fastq = f"/opt/conda2/pkgs/vsearch-2.4.3-0/bin/vsearch --fastq_filter {merged_outfile} " \
        f"--fastaout  {fasta} --fasta_width 0 --notrunclabels"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_fastq_fasta, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("##SBATCH -w, --nodelist=node01\n")
        handle.write("#SBATCH --mem=1000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{convert_fastq}\n")
        handle.write(f"echo {convert_unique_id}")
    os.chmod(str(run_fastq_fasta), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# Converting fastq to fasta:\n{convert_fastq}\n")

    cmd_fastq_fasta = f"sbatch -J {fastq_fasta_job_name} {run_fastq_fasta} --wait"
    try:
        fastq_fasta_slurm_id = subprocess.check_output(cmd_fastq_fasta, shell=True) \
            .decode(sys.stdout.encoding).strip()
        fastq_fasta_slurm_id = fastq_fasta_slurm_id.split(" ")[-1]
        fastq_fasta_slurm_out_file = pathlib.Path(chain_path, f"slurm-{fastq_fasta_slurm_id}.out")

    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"# vsearch fastq to fasta failed\n{e}\n")
        raise

    return fastq_fasta_slurm_out_file, convert_unique_id


def merged_files_gz(chain_path, script_folder, merged_outfile, itern, logfile):
    """

    :param chain_path:
    :param script_folder:
    :param merged_outfile:
    :param itern:
    :param logfile:
    :return:
    """
    # set job id
    gz_unique_id = uuid.uuid4()
    gzip_job_name = f'gzip{str(itern)}'
    run_gzip = pathlib.Path(script_folder, "run_gzip_merged.sh")
    cmd_gzip = f"gzip {merged_outfile}"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_gzip, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("#SBATCH -J gzip\n")
        handle.write("#SBATCH --mem=1000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{cmd_gzip}\n")
        handle.write(f"echo {gz_unique_id}")
    os.chmod(str(run_gzip), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# gzip on merged file:\n{cmd_gzip}\n")

    cmd_gzip = f"sbatch -J {gzip_job_name} {run_gzip} --parsable --wait"
    try:
        gz_merge_slurm_id = subprocess.check_output(cmd_gzip, shell=True) \
            .decode(sys.stdout.encoding).strip()
        gz_merge_slurm_id = gz_merge_slurm_id.split(" ")[-1]
        gz_merged_slurm_out_file = pathlib.Path(chain_path, f"slurm-{gz_merge_slurm_id}.out")

        zip_merged_file = pathlib.Path(f"{str(merged_outfile)}.gz")
        if not zip_merged_file.is_file():
            with open(logfile, "a") as handle:
                handle.write(f"# could not zip merged file\n# trying next sample\n")
            raise IOError

    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"# gzip on merged file failed\n{e}\n")
        raise e

    return gz_merged_slurm_out_file, gz_unique_id


def dereplicate_fasta(chain_path, scripts_folder, name_stem, derep_folder, fasta_to_derep_name_stem,
                      file_to_dereplicate, logfile):
    """

    :param chain_path:
    :param scripts_folder:
    :param name_stem:
    :param derep_folder:
    :param fasta_to_derep_name_stem:
    :param file_to_dereplicate:
    :param logfile:
    :return:
    """
    id_prefix = name_stem[3:6]
    unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
    # set job id
    derep_unique_id = uuid.uuid4()
    derep_job_name = f"{id_prefix}{unique_suffix}De"
    # set dereplicated outfile name
    dereplicated_file = pathlib.Path(derep_folder, f"{fasta_to_derep_name_stem}_unique.fasta")
    # create fastq to fasta SLURM file
    run_derep = pathlib.Path(scripts_folder, "run_derep.sh")
    derep_cmd = f"/opt/conda2/pkgs/vsearch-2.4.3-0/bin/vsearch --sizeout --derep_fulllength {file_to_dereplicate}" \
        f" --output {dereplicated_file} --fasta_width 0 --notrunclabels --threads 4"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_derep, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("##SBATCH -w, --nodelist=node01\n")
        handle.write("#SBATCH --mem=1000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{derep_cmd}\n")
        handle.write(f"echo {derep_unique_id}")
    os.chmod(str(run_derep), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# running dereplication\n{str(derep_cmd)}\n")

    cmd_derep = f"sbatch -J {derep_job_name} {run_derep} --ntasks=1 --cpus-per-task=4 --parsable --wait"
    try:
        derep_slurm_id = subprocess.check_output(cmd_derep, shell=True).decode(sys.stdout.encoding).strip()
        derep_slurm_id = derep_slurm_id.split(" ")[-1]
        derep_slurm_out_file = pathlib.Path(chain_path, f"slurm-{derep_slurm_id}.out")

        # fix file permissions for output
        os.chmod(str(dereplicated_file), 0o666)
    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"# vsearch dereplication failed\n{e}\n")
        raise

    return derep_slurm_out_file, derep_unique_id


def concat_fasta(chain_path, search_fasta_folder, scripts_folder, concated_outfile, logfile):
    # concatenate multiple fasta files
    cat_str = f"cat "
    for file in search_fasta_folder:
        cat_str += f"{str(file)} "
    # set job id
    cat_unique_id = uuid.uuid4()
    cat_job_name = f'concat'
    run_cat = pathlib.Path(scripts_folder, "run_cat.sh")
    concat_cmd = f"{cat_str} >> {concated_outfile}"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_cat, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("#SBATCH -J gzip\n")
        handle.write("#SBATCH --mem=1000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{concat_cmd}\n")
        handle.write(f"echo {cat_unique_id}")
    os.chmod(str(run_cat), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# concatenating multiple merged fasta files\n{concat_cmd}\n")

    cmd_cat = f"sbatch -J {cat_job_name} {run_cat} --parsable --wait"
    try:
        cat_slurm_id = subprocess.check_output(cmd_cat, shell=True).decode(sys.stdout.encoding).strip()
        cat_slurm_id = cat_slurm_id.split(" ")[-1]
        cat_slurm_out_file = pathlib.Path(chain_path, f"slurm-{cat_slurm_id}.out")

    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"# concatenating fasta files failed\n{e}\n")
        raise

    return cat_slurm_out_file, cat_unique_id


def check_slurm_jobs(job_name, job_list_to_check, sleep_time_sec, logfile):
    """

    :param job_name:
    :param job_list_to_check:
    :param sleep_time_sec:
    :param logfile:
    :return:
    """

    while True:
        if not job_list_to_check:
            break
        else:
            for indx, [job_file, unique_id] in enumerate(job_list_to_check):
                if job_file.is_file():
                    print(job_file)
                    input("enter")
                    with open(job_file, 'r') as slurm_out:
                        slurm_check = list(slurm_out)[-1].strip()
                    print(slurm_check)
                    input("enter")
                    if slurm_check == str(unique_id):
                        # remove the job entry from the list
                        del job_list_to_check[indx]
                        with open(logfile, "a") as handle:
                            handle.write(f"# {job_name}: job on {job_file} has finished\n")

            if job_list_to_check:
                time.sleep(sleep_time_sec)
    return True


def sonar_p1_call(chain_path, project_folder, scripts_folder, sample_name, sonar_version, dir_with_sonar1_files,
                  logfile, uniques_fasta):
    id_prefix = sample_name[3:6]
    unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
    # set job id
    sonar_p1_unique_id = uuid.uuid4()
    sonar1_job_name = f"{id_prefix}{unique_suffix}S1"
    # get sonar P1 version specific commands
    sonar_settings = {"heavy": f"python2 /opt/conda2/pkgs/sonar/annotate/1.1-blast_V.py "
                               f"-locus H -fasta {str(uniques_fasta)} "
                               f"-lib '/opt/conda2/pkgs/sonar/germDB/IgHJ.fa  "
                               f"-dlib /opt/conda2/pkgs/sonar/germDB/IgHD.fa  "
                               f"-clib /opt/conda2/pkgs/sonar/germDB/IgHC_CH1.fa -callFinal' -f",

                      "kappa": f"python /opt/conda2/pkgs/sonar/annotate/1.1-blast_V.py "
                               f"-locus K -fasta {str(uniques_fasta)}  "
                               f"-lib '/opt/conda2/pkgs/sonar/germDB/IgKJ.fa -noD -noC -callFinal' -f",

                      "lambda": f"python /opt/conda2/pkgs/sonar/annotate/1.1-blast_V.py "
                                f"-locus L -fasta {str(uniques_fasta)}  "
                                f"-lib '/opt/conda2/pkgs/sonar/germDB/IgLJ.fa -noD -noC -callFinal' -f"}

    sonar_version_setting = sonar_settings[sonar_version]
    run_sonar_p1 = pathlib.Path(scripts_folder, f"run_sonar_P1_{sonar_version}.sh")
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_sonar_p1, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("#SBATCH -w, --nodelist=bio-linux\n")
        handle.write("#SBATCH --mem=4000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{sonar_version_setting}\n")
        handle.write(f"echo {sonar_p1_unique_id}\n")
    os.chmod(str(run_sonar_p1), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# running Sonar P1 command from file:\n{str(sonar_version_setting)}\n")
    sonar1_run_cmd = f"sbatch -J {sonar1_job_name} {run_sonar_p1} --parsable"
    try:
        # change into sonar P1 targer directory
        os.chdir(dir_with_sonar1_files)
        sonar_p1_slurm_id = subprocess.check_output(sonar1_run_cmd, shell=True) \
            .decode(sys.stdout.encoding).strip()
        sonar_p1_slurm_id = sonar_p1_slurm_id.split(" ")[-1]
        sonar_p1_slurm_out_file = pathlib.Path(chain_path, f"slurm-{sonar_p1_slurm_id}.out")

    except subprocess.CalledProcessError as e:
        print("Sonar P1 encountered an error\ntrying next sample")
        with open(logfile, "a") as handle:
            handle.write(f"# Sonar P1 encountered \n{e}\n")
        os.chdir(project_folder)
        raise

    return sonar_p1_slurm_out_file, sonar_p1_unique_id


def sonar_p2_copy_files(chain_path, scripts_folder, dir_with_sonar1_output, dir_with_sonar1_work, target_folder,
                        logfile):
    """

    :param chain_path:
    :param scripts_folder:
    :param dir_with_sonar1_output:
    :param dir_with_sonar1_work:
    :param target_folder:
    :param logfile:
    :return:
    """
    # set job id
    sonar1_cp_unique_id = uuid.uuid4()
    sonar_p1_cp = f'snr1cp'
    run_snr1cp = pathlib.Path(scripts_folder, "run_sonar1_cp.sh")
    # copy files to desired directory
    cmd_copy_work = f"cp {dir_with_sonar1_work} {target_folder}"
    cmd_copy_output = f"cp {dir_with_sonar1_output} {target_folder}"
    slurm_outfile = str(pathlib.Path(chain_path, "slurm-%j.out"))
    with open(run_snr1cp, "w") as handle:
        handle.write("#!/bin/sh\n")
        handle.write("#SBATCH -J gzip\n")
        handle.write("#SBATCH --mem=1000\n")
        handle.write(f"#SBATCH -o {slurm_outfile}\n\n")
        handle.write(f"{cmd_copy_work}\n")
        handle.write(f"{cmd_copy_output}\n")
        handle.write(f"echo {sonar1_cp_unique_id}")
    os.chmod(str(run_snr1cp), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"# copyting Sonar P1 output to target directory\n{cmd_copy_work}"
                     f"\n{cmd_copy_output}")

    cmd_snr1_cp = f"sbatch -J {sonar_p1_cp} {run_snr1cp} --wait"
    try:
        sonar1_cp_slurm_id = subprocess.check_output(cmd_snr1_cp, shell=True) \
            .decode(sys.stdout.encoding).strip()
        sonar1_cp_slurm_id = sonar1_cp_slurm_id.split(" ")[-1]
        sonar1_cp_outfile = pathlib.Path(chain_path, f"slurm-{sonar1_cp_slurm_id}.out")

    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"Error copying Sonar P1 output to Sonar P2 folder\n{e}\n")
        raise

    return sonar1_cp_outfile, sonar1_cp_unique_id


def sonar_p2_call(chain_folder, sample_name, run_sonar2_trunc, known_mab_name, mab, scripts_folder, mab_name_file,
                  target_folder, parent_dir, logfile):
    """

    :param chain_folder:
    :param sample_name:
    :param run_sonar2_trunc:
    :param known_mab_name:
    :param mab:
    :param scripts_folder:
    :param mab_name_file:
    :param target_folder:
    :param parent_dir:
    :param logfile:
    :return:
    """
    # make Sonar P2 script
    id_prefix = sample_name[3:6]
    unique_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=2)).lower()
    # set job id
    sonar2_unique_id = uuid.uuid4()
    sonar2_job_name = f"{id_prefix}{unique_suffix}S2"
    if run_sonar2_trunc:
        sonar_p2_run = f"{known_mab_name}_{mab}_sonar_p2_run_trunc.sh"
        sonar_p2_run = pathlib.Path(scripts_folder, sonar_p2_run)
        sonar2_cmd = f"perl /opt/conda2/pkgs/sonar/lineage/2.1-calculate_id-div.pl -a {mab_name_file} " \
                     f"-g /opt/conda2/pkgs/sonar/germDB/IgHKLV_cysTruncated.fa -ap muscle -t 4"
        with open(sonar_p2_run, 'w') as handle:
            handle.write("#!/bin/sh\n")
            handle.write("##SBATCH -w, --nodelist=bio-linux\n")
            handle.write("#SBATCH --mem=4000\n")
            handle.write(f"#SBATCH -o {chain_folder}\n\n")
            handle.write(f"{sonar2_cmd}\n")
            handle.write(f"echo {sonar2_unique_id}\n")
    else:
        sonar_p2_run = f"{known_mab_name}_{mab}_sonar_p2_run.sh"
        sonar_p2_run = pathlib.Path(scripts_folder, sonar_p2_run)
        sonar2_cmd = f"perl /opt/conda2/pkgs/sonar/lineage/2.1-calculate_id-div.pl -a {mab_name_file} " \
                     f"-ap muscle -t 4"
        with open(sonar_p2_run, 'w') as handle:
            handle.write("#!/bin/sh\n")
            handle.write("##SBATCH -w, --nodelist=bio-linux\n")
            handle.write("#SBATCH --mem=4000\n\n")
            handle.write(f"{sonar2_cmd}\n")
            handle.write(f"echo {sonar2_unique_id}\n")
    os.chmod(str(sonar_p2_run), 0o777)

    with open(logfile, "a") as handle:
        handle.write(f"running sonar P2\n{sonar2_cmd}\n")

    sonar_p2_sumbit = f"sbatch -J {sonar2_job_name} {sonar_p2_run} --wait"
    try:
        # run sonar2
        os.chdir(str(target_folder))
        sonar2_cp_slurm_id = subprocess.check_output(sonar_p2_sumbit, shell=True) \
            .decode(sys.stdout.encoding).strip()
        sonar2_cp_slurm_id = sonar2_cp_slurm_id.split(" ")[-1]
        sonar2_outfile = pathlib.Path(parent_dir, f"slurm-{sonar2_cp_slurm_id}.out")
        os.chdir(parent_dir)
    except subprocess.CalledProcessError as e:
        with open(logfile, "a") as handle:
            handle.write(f"There was an error running sonar P2\n{e}\n")
        raise

    return sonar2_outfile, sonar2_unique_id


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
    sleep_time_sec = 60 * 2

    # check that raw files are gzipped
    check_gz_raw_jobs = []
    with open(logfile, "a") as handle:
        handle.write(f"# checking for raw files\n")
    for item in command_call_processing:
        sample_name = item[0]
        dir_with_raw_files = item[1]
        chain_path = dir_with_raw_files.parent
        script_folder = pathlib.Path(chain_path, "scripts")
        search_raw_files = list(dir_with_raw_files.glob(f"{sample_name}_R1.fastq.gz"))
        if not search_raw_files:
            with open(logfile, "a") as handle:
                handle.write(f"# No fastq.gz files in target folder for {sample_name}\n")

        search_raw_files_fastq = list(dir_with_raw_files.glob(f"{sample_name}_R1.fastq"))
        if search_raw_files_fastq:
            with open(logfile, "a") as handle:
                handle.write(f"# unzipped files sound\n# gziping the files\n")
            gz_slurm_out_file, gz_unique_id = raw_files_gz(chain_path, sample_name, dir_with_raw_files,
                                                           script_folder, logfile)
            check_gz_raw_jobs.append([gz_slurm_out_file, gz_unique_id])

    # check for completion of gzip submissions
    if check_gz_raw_jobs:
        print("waiting for gzip on raw files")
        with open(logfile, "a") as handle:
            handle.write(f"# waiting for gzip on raw files\n")
        check_slurm_jobs("gzip of raw files", check_gz_raw_jobs, sleep_time_sec, logfile)

    # run pear on raw files
    print("running PEAR")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running PEAR\n\n")
    check_pear_jobs = []
    merged_folders = []
    for item in command_call_processing:
        sample_name = item[0]
        with open(logfile, "a") as handle:
            handle.write(f"\n# Processing {sample_name}\n")
        print(f"\n# Processing {sample_name}\n")
        dir_with_raw_files = item[1]
        chain_path = dir_with_raw_files.parent
        script_folder = pathlib.Path(chain_path, "scripts")
        merged_folder = pathlib.Path(chain_path, "2_merged_filtered")
        fasta_folder = pathlib.Path(chain_path, "3_fasta")
        derep_folder = pathlib.Path(chain_path, "4_dereplicated")

        # remove old merged files if present
        merged_search = list(merged_folder.glob(f"{sample_name}*.fastq*"))
        if merged_search:
            with open(logfile, "a") as handle:
                handle.write(f"# removing old merged file for {sample_name}\n")
            for file in merged_search:
                os.unlink(str(file))

        # remove old fasta file
        with open(logfile, "a") as handle:
            handle.write(f"removing old fasta file for {sample_name}\n")
        fasta_search = list(fasta_folder.glob(f"{sample_name}*.fasta*"))
        if fasta_search:
            for file in fasta_search:
                os.unlink(str(file))

        # remove old dereplicated file
        with open(logfile, "a") as handle:
            handle.write(f"removing old dereplicated file for {sample_name}\n")
        derep_search = list(derep_folder.glob("*_unique.fasta"))
        if derep_search:
            for file in derep_search:
                os.unlink(str(file))

        # check for zipped files
        search_zip_raw_files = list(dir_with_raw_files.glob(f"{sample_name}_R*.fastq.gz"))
        if search_zip_raw_files:
            print("raw files found")
            with open(logfile, "a") as handle:
                handle.write(f"# raw files found in target folder {str(search_zip_raw_files)}\n")

            for i, file_r1 in enumerate(search_zip_raw_files):
                file_r2 = str(file_r1).replace("_R1.fastq.gz", "_R2.fastq.gz")
                file_r2 = pathlib.Path(file_r2)
                # if can't find the matched paired R2 file, skip and go to next R1 file
                if not file_r2.is_file():
                    print(f"R2 paired file not found\nSkipping this sample: {file_r1}")
                    with open(logfile, "a") as handle:
                        handle.write(f"# No R2 paired file not found\n# Skipping this sample: {file_r1}\n")
                    continue
                # if both R1 and R2 files are present, continue with the processing
                else:
                    # get for free disk space
                    free_space, percent_free = disk_space_checker(path)
                    file_size = os.path.getsize(file_r1) / gb * 2
                    if free_space < file_size * 4:
                        print(f"not enough disk space to process the raw files\nNeed{file_size * 4}Gb\n"
                              f"you have{free_space}Gb free")
                        with open(logfile, "a") as handle:
                            handle.write(
                                f"# not enough disk space to process the raw files\n# Need {file_size * 4}Gb\n"
                                f"# you have{free_space}Gb free\n")
                        sys.exit("exiting")

                    merged_file_name = pathlib.Path(merged_folder, sample_name)
                    pear_slurm_out_file, pear_unique_id = \
                        run_pear(chain_path, sample_name, script_folder, file_r1, file_r2, merged_file_name, i, logfile)

                    # collect submission
                    check_pear_jobs.append([pear_slurm_out_file, pear_unique_id])
                    merged_folders.append(merged_folder)

    # check for completion of pear submissions
    print("waiting for pear")
    with open(logfile, "a") as handle:
        handle.write(f"# waiting for pear merging of reads\n")
    pear_sleep = sleep_time_sec * 10
    check_slurm_jobs("pear merge", check_pear_jobs, pear_sleep, logfile)

    # remove unmerged files
    for merged_folder in merged_folders:
        unassmebled_search = pathlib.Path(merged_folder).glob(f"{sample_name}*unassembled*.fastq")
        discarded_search = pathlib.Path(merged_folder).glob(f"{sample_name}*discarded.fastq")
        for file in unassmebled_search:
            os.unlink(str(file))
        for file in discarded_search:
            os.unlink(str(file))

    # convert to fasta
    print("running fastq 2 fasta")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running fastq 2 fasta\n\n")
    check_convert_jobs = []
    fasta_file_list = []
    for j, item in enumerate(command_call_processing):
        sample_name = item[0]
        print(f"# Processing {sample_name}\n")
        with open(logfile, "a") as handle:
            handle.write(f"# Processing {sample_name}\n")
        dir_with_raw_files = item[1]
        chain_path = dir_with_raw_files.parent
        script_folder = pathlib.Path(chain_path, "scripts")
        merged_folder = pathlib.Path(chain_path, "2_merged_filtered")
        fasta_folder = pathlib.Path(chain_path, "3_fasta")
        meged_outfile = pathlib.Path(merged_folder, f"{sample_name}.assembled.fastq")
        fasta = f"{str(meged_outfile.stem)}.fasta"
        fasta = pathlib.Path(fasta_folder, fasta)

        # fix file permissions for output
        os.chmod(str(meged_outfile), 0o666)
        fastq_fasta_slurm_out_file, convert_unique_id = \
            fastq2fasta(chain_path, sample_name, script_folder, fasta, meged_outfile, j, logfile)

        # collect submission
        check_convert_jobs.append([fastq_fasta_slurm_out_file, convert_unique_id])
        fasta_file_list.append(fasta)

    # check for completion of pear submissions
    print("waiting for fastq to fasta conversion")
    with open(logfile, "a") as handle:
        handle.write(f"# waiting for fastq to fasta conversion\n")
    check_slurm_jobs("fastq2fasta", check_convert_jobs, sleep_time_sec, logfile)

    # gzip the merged files
    print("running gzip on merged files")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running gzip on merged files\n\n")
    check_gz_merged_jobs = []
    for j, item in enumerate(command_call_processing):
        sample_name = item[0]
        print(f"# Processing {sample_name}")
        with open(logfile, "a") as handle:
            handle.write(f"# Processing {sample_name}\n")
        dir_with_raw_files = item[1]
        chain_path = dir_with_raw_files.parent
        script_folder = pathlib.Path(chain_path, "scripts")
        merged_folder = pathlib.Path(chain_path, "2_merged_filtered")
        fasta_folder = pathlib.Path(chain_path, "3_fasta")
        meged_outfile = pathlib.Path(merged_folder, f"{sample_name}.assembled.fastq")
        fasta = f"{str(meged_outfile.stem)}.fasta"
        fasta = pathlib.Path(fasta_folder, fasta)
        # compress pear if the merged file was successfully converted to a fasta file
        if fasta.is_file():
            # fix file permissions for output
            os.chmod(str(fasta), 0o666)
            print("gzip on merged fastq file")
            with open(logfile, "a") as handle:
                handle.write(f"# gzipping merged fastq file\n")
            gz_merged_slurm_out_file, gz_unique_id = merged_files_gz(chain_path, script_folder, meged_outfile, j,
                                                                     logfile)
            # collect submission
            check_gz_merged_jobs.append([gz_merged_slurm_out_file, gz_unique_id])

    # check for completion of pear submissions
    print("waiting for gzip on merged files")
    with open(logfile, "a") as handle:
        handle.write(f"# waiting for gzip on merged files\n")
    check_slurm_jobs("gzip on merged", check_gz_merged_jobs, sleep_time_sec, logfile)

    # collect all the files that will be dereplicated
    files_to_derep_dict = collections.defaultdict(list)
    for item in command_call_processing:
        sample_name = item[0]
        dir_with_raw_files = item[1]
        chain_path = dir_with_raw_files.parent
        files_to_derep_dict[str(chain_path)].append(sample_name)

    # concat (if needed) and dereplicate fasta files
    print("running dereplication (and concatenation if needed)")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running dereplication (and concatenation if needed)\n\n")
    concat_files = []
    check_derep_jobs = []
    for chain_folder, sample_name_list in files_to_derep_dict.items():
        check_concat_jobs = []
        scripts_folder = pathlib.Path(chain_folder, "scripts")
        fasta_folder = pathlib.Path(chain_folder, "3_fasta")
        derep_folder = pathlib.Path(chain_folder, "4_dereplicated")
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
            print("concatenating multiple merged.fasta files")
            with open(logfile, "a") as handle:
                handle.write(f"concatenating multiple merged.fasta files\n")
            concated_outfile = pathlib.Path(fasta_folder, f"{fasta_to_derep_name_stem}_concatenated.fasta")
            if concated_outfile.is_file():
                # remove existing file so that you don't accidentally concatenate in duplicate
                print(f"{concated_outfile}\nalready exists\nthis file will be overwritten")
                with open(logfile, "a") as handle:
                    handle.write(f"{concated_outfile}\nalready exists\nthis file will be overwritten/n")
                os.unlink(str(concated_outfile))

            # concatenate multiple fasta files
            cat_slurm_out_file, cat_unique_id = concat_fasta(chain_folder, search_fasta_folder, scripts_folder,
                                                             concated_outfile, logfile)
            check_concat_jobs.append([cat_slurm_out_file, cat_unique_id])

            # check for completion of concat job
            check_slurm_jobs("concat fasta", check_concat_jobs, sleep_time_sec, logfile)

            file_to_dereplicate = concated_outfile
        # only one fasta file for this sample
        elif len(search_fasta_folder) == 1:
            # copy fasta to 4_dereplicated
            file_to_dereplicate = search_fasta_folder[0]
        else:
            print(f"no Fasta files were found\nskipping this sample: {sample_name_list}")
            with open(logfile, "a") as handle:
                handle.write(f"# no Fasta files were found\nskipping this sample: {sample_name_list}")
            continue

        # dereplicate sequences using vsearch
        print(f"dereplicating fasta file {file_to_dereplicate}")
        with open(logfile, "a") as handle:
            handle.write(f"\n# dereplicating file {file_to_dereplicate}\n")

        derep_slurm_out_file, derep_unique_id = dereplicate_fasta(chain_folder, scripts_folder, name_stem,
                                                                  derep_folder, fasta_to_derep_name_stem,
                                                                  file_to_dereplicate, logfile)
        check_derep_jobs.append([derep_slurm_out_file, derep_unique_id])
        concat_files.append(fasta_folder)

    # check for derep completion
    print("waiting for dereplication on files")
    with open(logfile, "a") as handle:
        handle.write(f"# waiting for dereplication on files\n")
    check_slurm_jobs("dereplication of fasta", check_derep_jobs, sleep_time_sec, logfile)

    # remove the concatenated files
    print("removing any concatenated files")
    with open(logfile, "a") as handle:
        handle.write(f"# removing any concatenated files\n")
    for fasta_folder in concat_files:
        # remove non-dereplicated file if you had to concatenate multiple files
        for file in pathlib.Path(fasta_folder).glob(f"*_concatenated.fasta"):
            print(f"removing concatenated fasta file: {file}")
            with open(logfile, "a") as handle:
                handle.write(f"# removing {file}\n")
            os.unlink(str(file))

    with open(logfile, "a") as handle:
        handle.write(f"# processing of sample completed\n")


def step_2_run_sonar_p1(command_call_sonar_1, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_1: list of samples to run sonar1 on
    :param logfile: (str) path and name of the logfile
    :return:
    """
    gb = (1024 * 1024) * 1024
    slurm_submission_jobs = []
    sleep_time_sec = 60 * 10 * 6

    print("running sonar P1")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running sonar P1\n")
    for item in command_call_sonar_1:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        chain_folder = pathlib.Path(dir_with_sonar1_files.parent)
        scripts_folder = pathlib.Path(dir_with_sonar1_files.parent, "scripts")
        dir_with_sonar1_work = pathlib.Path(dir_with_sonar1_files, "work")
        dir_with_sonar1_output = pathlib.Path(dir_with_sonar1_files, "output")

        with open(logfile, "a") as handle:
            handle.write(f"\n{'-'*10}\n# running Sonar P1 on {sample_name}\n")

        # check whether there is already sonar P1 output in the target folders (if this is a rerun)
        search_work = list(dir_with_sonar1_work.glob("*.*"))
        search_output = list(dir_with_sonar1_output.glob("*.*"))
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
        os.chdir(project_folder)

        # check that only one file in target dir
        search_derep_fastas = list(dir_with_sonar1_files.glob("*.fasta"))
        if len(search_derep_fastas) > 1:
            print(f"multiple fasta files found in {dir_with_sonar1_files}\n"
                  f"Sonar1 can only run if there is one fasta file in this folder\ntrying next sample")
            with open(logfile, "a") as handle:
                handle.write(f"# multiple fasta files in {dir_with_sonar1_files}\n# trying next sample")
            continue
        elif len(search_derep_fastas) == 0:
            print(f"No fasta files found in {dir_with_sonar1_files}\n"
                  f"Sonar1 can only run if there is one fasta file in this folder\ntrying next sample")
            with open(logfile, "a") as handle:
                handle.write(f"# No fasta files in {dir_with_sonar1_files}\n")
        else:
            # check that you have enough space
            uniques_fasta = search_derep_fastas[0]

            # fix file permissions for output
            os.chmod(str(uniques_fasta), 0o666)
            # get for free disk space
            free_space, percent_free = disk_space_checker(project_folder)
            fasta_size = os.path.getsize(uniques_fasta) / gb
            if free_space < fasta_size * 10:
                print(f"not enough free space to run sonar P1\nFree space is {free_space}Gb"
                      f"# you expected to need up to {fasta_size * 10}Gb for sonar P1")
                with open(logfile, "a") as handle:
                    handle.write(
                        f"# Not enough free space to sonar P1\nFree space is {free_space}Gb {percent_free}%\n"
                        f"# you expected to need up to {fasta_size * 10}Gb for sonar P1")
                sys.exit("exiting")

            # make sonar p1 call
            sonar_p1_slurm_out_file, sonar_p1_unique_id = sonar_p1_call(chain_folder, project_folder, scripts_folder,
                                                                        sample_name, sonar_version,
                                                                        dir_with_sonar1_files, uniques_fasta, logfile)
            slurm_submission_jobs.append([sonar_p1_slurm_out_file, sonar_p1_unique_id])

    # search for finished sonar P1 jobs and remove from list, Holds pipeline until Sonar P1 jobs are done
    check_slurm_jobs("sonar P1", slurm_submission_jobs, sleep_time_sec, logfile)


def step_3_run_sonar_2(command_call_sonar_2, fasta_sequences, run_sonar2_trunc, logfile):
    """
    function to automate calling sonar1 on NGS Ab dereplicated fasta files
    :param command_call_sonar_2: list of samples to run sonar2 on
    :param fasta_sequences: (dict) dict of all the mAb fasta sequences (full and crd3)
    :param run_sonar2_trunc: (Bool) True if flag is set, default is False, will run sonar_trunc call if True
    :param fasta_sequences: (dict) dict containing key = seq name, value = sequence for all known mAb sequences
    :param logfile: (str) path and name of the logfile
    :return:
    """
    gb = (1024 * 1024) * 1024
    sleep_time_sec = 60 * 10 * 2

    print("running sonar P2")
    with open(logfile, "a") as handle:
        handle.write(f"\n# {'-' * 10}\n# running sonar P2\n")

    slurm_sonar2_submission_jobs = []
    for item in command_call_sonar_2:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        chain_folder = pathlib.Path(dir_with_sonar1_files.parent)
        scripts_folder = pathlib.Path(dir_with_sonar1_files.parent, "scripts")
        parent_dir = dir_with_sonar1_files.parents[3]
        chain = item[2]
        known_mab_name = item[3]
        target_folder_full_ab = pathlib.Path(f"5_{known_mab_name}", "full_ab")
        target_folder_crdh3 = pathlib.Path(f"5_{known_mab_name}", "crd3")
        fullab_name = known_mab_name + f"_{chain}_fullab"
        cdr3_name = known_mab_name + f"_{chain}_cdr3"
        try:
            mab_sequence_fullab = fasta_sequences[fullab_name]
            mab_sequence_cdr3 = fasta_sequences[cdr3_name]
        except KeyError:
            with open(logfile, "a") as handle:
                handle.write(f"# mab names in settings file don't match the names in the fasta file\n")
            raise

        # check for mab_sequences folder
        fasta_sequences = pathlib.Path(parent_dir, "mab_sequences")
        if not fasta_sequences.is_dir():
            print(f"'mab_sequences' folder not found\n"
                  f"You need to create a folder 'mab_sequences' in the project folder, ie:\n{fasta_sequences}"
                  f"\nand copy your fasta file in there")
            print("making the folder for you")
            with open(logfile, "a") as handle:
                handle.write(f"# mab_sequences' folder not found\n"
                             f"# You need to create a folder 'mab_sequences' in the project folder, ie:\n"
                             f"{fasta_sequences}\n# making the folder for you")
            fasta_sequences.mkdir(mode=0o777, parents=True, exist_ok=True)
            sys.exit("exiting")

        # check for Sonar P1 files
        dir_with_sonar1_work = pathlib.Path(dir_with_sonar1_files, "work")
        dir_with_sonar1_output = pathlib.Path(dir_with_sonar1_files, "output")
        for sonar1_dir in [dir_with_sonar1_work, dir_with_sonar1_output]:
            sonar_p1_search = sonar1_dir.glob("*")
            if not sonar_p1_search:
                print(f"No Sonar P1 files detected in: {dir_with_sonar1_files}")
                with open(logfile, "a") as handle:
                    handle.write(f"# No Sonar P1 files detected in: {dir_with_sonar1_files}\n")
                    sys.exit()

        # run on full antibody and on cdr3 region
        to_run = [(target_folder_full_ab, fullab_name, mab_sequence_fullab),
                  (target_folder_crdh3, cdr3_name, mab_sequence_cdr3)]
        mab = ["fullab", "cdr3"]

        for i, (target_folder, mab_name, mab_seq) in enumerate(to_run):
            slurm_cp_submission_jobs = []
            # make fasta file for ab run
            mab_name_file = f"{item}.fasta"
            mab_name_file = pathlib.Path(target_folder, mab_name_file)
            with open(mab_name_file, 'w') as handle:
                handle.write(f">{mab_name}\n{mab_sequence_fullab}\n")
            os.chmod(str(mab_name_file), 0o666)
            # set target folders
            dir_with_sonar2_work = pathlib.Path(target_folder, "work")
            dir_with_sonar2_output = pathlib.Path(target_folder, "output")
            # check that target dirs have unmodified sonar P1 output
            sonar_p1_alread_cp_work = is_same(dir_with_sonar1_work, dir_with_sonar2_work)
            sonar_p1_alread_cp_output = is_same(dir_with_sonar1_output, dir_with_sonar2_output)
            copy_sonar_p1 = True
            for sonar2_dir in [dir_with_sonar2_work, dir_with_sonar2_output]:
                sonar_p2_search = sonar2_dir.glob(".*")
                # if sonar P2 folders not empty
                if sonar_p2_search:
                    # if they are identical to sonar P1 folders don't allow copy
                    if sonar_p1_alread_cp_work and sonar_p1_alread_cp_output:
                        print("unmodified sonar P1 output found in sonar P2 target folder\nNot re-copying the files")
                        with open(logfile, "a") as handle:
                            handle.write(f"# unmodified sonar P1 output found in sonar P2 target folder\n"
                                         f"# Not re-copying the files\n")
                        copy_sonar_p1 = False
                    # if not identical, remove files and allow copy
                    else:
                        print(f"Modified sonar P1 output detected in sonar P2 target folder {item}\ndeleting files")
                        with open(logfile, "a") as handle:
                            handle.write(f"# Modified sonar P1 output detected in sonar P2 target folder {item}\n"
                                         f"# deleting files\n")
                        for file in sonar_p2_search:
                            os.unlink(str(file))

                # copy data if sonar P1 output not already (unmodified) in target dir
                if copy_sonar_p1:
                    # check for free disk space
                    free_space, percent_free = disk_space_checker(parent_dir)
                    sonar_p1_size = round(folder_size_checker(str(dir_with_sonar1_output)) / gb, 3)
                    if free_space < sonar_p1_size * 3:
                        print(f"not enough free space to run sonar P2\nFree space is {free_space}Gb"
                              f"# you expected to need up to {sonar_p1_size * 3}Gb for sonar P2")
                        with open(logfile, "a") as handle:
                            handle.write(f"# Not enough free space to sonar P2\n"
                                         f"# Free space is {free_space}Gb {percent_free}%\n"
                                         f"# you expected to need up to {sonar_p1_size * 3}Gb for sonar P2")
                        sys.exit("exiting")

                    sonar1_cp_outfile, sonar1_cp_unique_id = sonar_p2_copy_files(chain_folder, scripts_folder,
                                                                                 dir_with_sonar1_output,
                                                                                 dir_with_sonar1_work,
                                                                                 target_folder, logfile)
                    # collect cp job details
                    slurm_cp_submission_jobs.append([sonar1_cp_outfile, sonar1_cp_unique_id])

                    # check slurm sonar p1 copy jobs
                    check_slurm_jobs("copy sonar P1 files", slurm_cp_submission_jobs, sleep_time_sec, logfile)

                    # check the files copied correctly
                    sonar_p1_alread_cp_work = is_same(dir_with_sonar1_work, dir_with_sonar2_work)
                    sonar_p1_alread_cp_output = is_same(dir_with_sonar1_output, dir_with_sonar2_output)
                    if not sonar_p1_alread_cp_output and sonar_p1_alread_cp_work:
                        print("Copied files are different from sonar P1 original files\nCopy must have failed")
                        with open(logfile, "a") as handle:
                            handle.write(f"# Copied files are different from sonar P1 original files\n"
                                         f"# Copy must have failed")
                        sys.exit("exiting")

            sonar2_outfile, sonar2_unique_id = sonar_p2_call(chain_folder, sample_name, run_sonar2_trunc,
                                                             known_mab_name, mab[i], scripts_folder, mab_name_file,
                                                             target_folder, parent_dir, logfile)

            slurm_sonar2_submission_jobs.append([sonar2_outfile, sonar2_unique_id])

    # search for finished sonar P2 jobs and remove from list, Holds pipeline until Sonar P2 jobs are done
    check_slurm_jobs("sonar P2", slurm_sonar2_submission_jobs, sleep_time_sec, logfile)


def main(path, settings, fasta_file=None, run_sonar2_trunc=False):
    """
    create folder structure for a sequencing project
    :param path: (str) path to where the folders should be made
    :param settings: (str) the path and csv file with the run settings
    :param fasta_file: (str) the path and name of the fasta file with your mAb sequences
    :param run_sonar2_trunc: (Bool) set flag to true to run sonar_trunc
    """

    path = pathlib.Path(path).absolute()
    if not path.is_dir():
        print("path to project folder is not a directory\n")
        sys.exit("exiting")

    # change the cwd to project path
    os.chdir(path)

    # check disk space
    free_space, percent_free = disk_space_checker(path)
    print(f"Free disk space remaining: {free_space}Gb ({percent_free}%)")

    # set paths and initialize the log file
    project_name = path.parts[-1]
    settings = pathlib.Path(settings).absolute()

    if not settings.is_file():
        print("settings file not found\n")
        sys.exit("exiting")

    settings_dataframe = pd.read_csv(settings, sep=None, engine='python')
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(path, f"{time_stamp}_log_file.txt")
    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for project: {project_name}\n")
        handle.write(f"# Free space at start = {free_space}Gb ({percent_free}%)\n\n")
    os.chmod(str(log_file), 0o777)

    # check that settings file has the correct headings
    settings_dict = settings_checker(settings_dataframe, log_file)

    # make a list of the jobs to run and the settings required and make the required folders
    list_all_settings = extract_settings_make_folders(path, settings_dict, log_file)

    if not list_all_settings:
        print("No jobs were found in the settings file\nCheck that the file format was correct "
              "and that it contains entries")
        with open(log_file, "w") as handle:
            handle.write(f"# No jobs were found in the settings file\n"
                         f"Check that the file format was correct and that it contains entries\n")
        sys.exit("exiting")

    # generate the job lists
    command_call_processing, command_call_sonar_1, command_call_sonar_2 = make_job_lists(path, list_all_settings,
                                                                                         log_file)
    # only run sample processing if one or more files were specified
    if command_call_processing:
        print("processing jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# processing jobs found\n")
        # move files from 0_new_data, into correct directory, if necessary
        move_raw_data(path, settings_dataframe, log_file)

        # get rid of duplicate entries in sample processing list, if present
        command_call_processing = list(set(tuple(x) for x in command_call_processing))
        step_1_run_sample_processing(path, command_call_processing, log_file)
    else:
        print("no processing jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# no processing jobs found\n")

    # only run sonar1 if one or more files were specified
    if command_call_sonar_1:
        print("sonar P1 jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# sonar P1 jobs found\n")
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
    else:
        print("no sonar P1 jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# no sonar P1 jobs found\n")
    # only run sonar 2 if one or more files were specified
    if command_call_sonar_2:
        print("sonar P2 jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# sonar P2 jobs found\n")
        if fasta_file:
            fasta_sequences = fasta_to_dct(fasta_file)
            step_3_run_sonar_2(command_call_sonar_2, fasta_sequences, run_sonar2_trunc, log_file)
        else:
            print("no fasta file specified for Sonar P2 to run\n")
            sys.exit("exiting")
    else:
        print("no sonar P2 jobs found")
        with open(log_file, "w") as handle:
            handle.write(f"# no sonar P2 jobs found\n")
    print("Done")


if __name__ == "__main__":
    args = docopt(__doc__, version='V0.1.0')

    main(path=args["<project_path>"], settings=args["<settings_file>"], fasta_file=args["<fasta_file>"],
         run_sonar2_trunc=args["--run_sonar2_trunc"])
