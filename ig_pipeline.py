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
import sys
import os
from docopt import docopt
import pandas as pd
import pathlib
import subprocess


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


def step_0_make_folders(path, lineage, time_point, chain, primer_name):
    """
    function to create the nested folder structure for the nAb pileline, if they don't already exist
    :param path: (str) the path to the project folder
    :param lineage: (str) Ab lineage name
    :param time_point: (str) the time point the NGS data was sampled from
    :param chain: (str) the Ab chain (heavy, kappa, lambda)
    :param primer_name: (str) the primer name used for target amplification
    :return: None
    """

    primer_name = "5_" + primer_name
    pathlib.Path(path, lineage, time_point, chain, "scripts").mkdir(mode=0o777, parents=True,
                                                                       exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "1_raw_data").mkdir(mode=0o777, parents=True,
                                                                       exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "2_merged_filtered").mkdir(mode=0o777, parents=True,
                                                                              exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "3_fasta").mkdir(mode=0o777, parents=True,
                                                                    exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "4_dereplicated", "work").mkdir(mode=0o777, parents=True,
                                                                                   exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, "4_dereplicated", "output").mkdir(mode=0o777, parents=True,
                                                                                     exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, primer_name, "crdh3", "work").mkdir(mode=0o777, parents=True,
                                                                                       exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, primer_name, "crdh3", "output").mkdir(mode=0o777, parents=True,
                                                                                         exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, primer_name, "full_ab", "work").mkdir(mode=0o777, parents=True,
                                                                                         exist_ok=True)
    pathlib.Path(path, lineage, time_point, chain, primer_name, "full_ab", "output").mkdir(mode=0o777, parents=True,
                                                                                           exist_ok=True)


def unzip_files(path):
    """
    a function to upzip any .gz or .zip compressed files
    :param path: the path to the 0_new_data folder
    :return:
    """
    search = pathlib.Path(path, "0_new_data").glob("*")

    for file in search:

        print(file)
        # do a '-' t0 '_' rename
        # cmd = f"rename -v '-' '_' {file}"
        cmd = f"rename 's/-/_/g' {file}"
        subprocess.call(cmd, shell=True)

        # uncompress the files
        suffix = file.suffix
        print(suffix)
        if suffix == ".zip":
            cmd = f"unzip -o -d ./ {file}"
            subprocess.call(cmd, shell=True)
            unzipped = str(file).replace(".zip", "")

            # remove the .zip file after unzipping
            if os.path.isfile(unzipped):
                cmd = f"rm {file}"
                subprocess.call(cmd, shell=True)
            else:
                print(f"unzipping of {file} failed")

        elif suffix == ".gz":
            cmd = f"gunzip {file}"
            subprocess.call(cmd, shell=True)

        elif suffix == ".fastq":
            print("the file is already unziped")

        else:
            print(f"unknown file ending for {file}\nExpected .fastq, .zip or .gz")


def move_raw_data(path, settings_dataframe):
    """
    function to move raw data to the correct folder within the project
    :param path: path to the project folder
    :param settings_dataframe: the dataframe object form the settings file
    :return: a list of all the files that were moved
    """
    raw_files = pathlib.Path(path, "0_new_data").glob("*.fastq")
    n = 0
    # rename raw files
    for n, file in enumerate(raw_files):
        new_name = str(file).replace("-", "_")
        os.rename(str(file), new_name)

    # find where to move each raw file and move it to the right folder
    raw_files = pathlib.Path(path, "0_new_data").glob("*.fastq")
    for n, file in enumerate(raw_files):
        full_name = file
        name = file.stem
        parts = name.split("_")
        search_name = "_".join(parts[:5])
        copy_name = "_".join(parts[:5]) + f"_{parts[7]}.fastq"
        fields = settings_dataframe.loc[settings_dataframe['sample_name'] == search_name].head(1)
        if fields.empty:
            print(f"Your sample name {search_name}\nwas not found in the settings 'sample_name' column\n"
                  f"Please fix the file name or settings file accordingly")
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
    if n == 0:
        print("No fastq (fastq.gz/fastq.zip) files in 0new_data\n")


def step_1_run_sample_processing(command_call_processing):
    for item in command_call_processing:
        sample_name = item[0]
        dir_with_raw_files = item[1]
        parent_path = dir_with_raw_files.parent

        merged_folder = pathlib.Path(parent_path, "2_merged_filtered")
        fasta_folder = pathlib.Path(parent_path, "3_fasta")
        derep_folder = pathlib.Path(parent_path, "4_dereplicated")

        search_raw_files = dir_with_raw_files.glob(f"{sample_name}*_R1_*.fastq*")
        for file_R1 in search_raw_files:
            file_R2 = str(file_R1).replace("_R1_", "_R2_")
            file_R2 = pathlib.Path(file_R2)
            if not file_R2.is_file():
                print(f"R2 paired file not found\nSkipping this sample: {file_R1}")
            else:
                # check if file is zipped, unzip if needed
                extension_R1 = file_R1.suffix
                extension_R2 = file_R2.suffix
                if extension_R1 == ".gz":
                    cmd_unzip_R1 = f"gunzip {file_R1}"
                    subprocess.call(cmd_unzip_R1, shell=True)
                    file_R1 = pathlib.Path(str(file_R1).replace(".gz", ""))
                if extension_R2 == ".gz":
                    cmd_unzip_R2 = f"gunzip {file_R2}"
                    subprocess.call(cmd_unzip_R2, shell=True)
                    file_R2 = pathlib.Path(str(file_R2).replace(".gz", ""))

                # run pear
                pear_outfile = file_R1.name.replace("fastq", "")
                pear_outfile = f"{pear_outfile}_merged.fastq"
                merged_file_name = pathlib.Path(merged_folder, pear_outfile)
                pear_job_name = ''
                run_pear = pathlib.Path(parent_path, "scripts", "run_pear.sh")
                with open(run_pear, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=node01\n")
                    handle.write("#SBATCH --mem=1000\n")
                    handle.write(f"pear -f {file_R1} -r {file_R2} -o {merged_file_name} -p 0.001 -j 4")

                cmd_pear = f"sbatch -J {pear_job_name} {run_pear}"
                subprocess.call(cmd_pear, shell=True)

                # remove unmerged files?
                unassmebled_search = pathlib.Path(merged_folder).glob(f"{sample_name}*unassembled*.fastq")
                discarded_search = pathlib.Path(merged_folder).glob(f"{sample_name}*discarded.fastq")
                for file in unassmebled_search:
                    cmd_rm = f"rm {file}"
                    subprocess.call(cmd_rm, shell=True)
                for file in discarded_search:
                    cmd_rm = f"rm {file}"
                    subprocess.call(cmd_rm, shell=True)

                # compress 1_raw_data
                cmd_zip = f"gzip {file_R1} {file_R2}"
                subprocess.call(cmd_zip, shell=True)

                # convert to fasta
                seqmagick_job_name = ''
                fastq = merged_file_name
                fasta = str(fastq.name).replace("fastq", "fasta")
                fasta = pathlib.Path(fasta_folder, fasta)
                run_seqmagick = pathlib.Path(parent_path, "scripts", "run_convert_to_fasta.sh")
                with open(run_seqmagick, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    handle.write("##SBATCH -w, --nodelist=node01\n")
                    handle.write("#SBATCH --mem=1000\n")
                    handle.write(f"seqmagick convert {fastq} {fasta}")

                cmd_seqmagick = f"sbatch -J {seqmagick_job_name} {run_seqmagick}"
                subprocess.call(cmd_seqmagick, shell=True)

                # compress pear
                cmd_zip = f"gzip {}"
                subprocess.call(cmd_zip, shell=True)

        # concatenate multiple files
        # check if have multiple files

        search_fasta_files = pathlib.Path(fasta_folder).glob("*.fasta")
        all_fasta = [x for x in search_fasta_files]
        cat_file = sample_name + '_concatenated.fastq'
        if len(all_fasta) > 1:
            for file in files_to_cat:
                concat_cmd = f"cat {file} {cat_file}"
        else:
            # copy fasta to 4_dereplicated
        # zip fasta files

        # dereplicate sequences using vsearch
        derep_file = ''
        dereplicate_cmd = f"vsearch dereplicate {cat_file} {derep_file}"

        # rename to sample_name_unique.fasta
        out_name = derep_file.replace(".fasta", "_unique.fasta")
        rename_cmd = f"rename -v {derep_file} {out_name}"
        subprocess.call(rename_cmd, shell=True)


def step_2_run_sonar_1(command_call_sonar_1):
    for item in command_call_sonar_1:
        sample_name = item[0]
        dir_with_sonar1_files = item[1]
        sonar_version = item[2]

        # check that only one file in target dir

        job_name = ''

        if sonar_version == "heavy":
            cmd = f"sbatch -J {job_name} /home/job_scripts/sonar/annotate/run_sonar-P1-H.sh"
        elif sonar_version == "kappa":
            cmd = f"sbatch -J {job_name} /home/job_scripts/sonar/annotate/run_sonar_light.sh"
        elif sonar_version == "lambda":
            cmd = f"sbatch -J {job_name} /home/job_scripts/sonar/annotate/run_sonar_P1-kappa.sh"

        # tar.gz output


def step_3_run_sonar_2(command_call_sonar_2):
    for item in command_call_sonar_2:
        sample_name = item[0]
        dir_with_sonar2_files = item[1]
        sonar_version = item[2]
        primer_name = item[3]
        known_mab_name = item[4]

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
    get_script_path = pathlib.Path(__file__)
    get_script_path = get_script_path.absolute()
    script_folder = get_script_path.parent
    print(script_folder)

    path = pathlib.Path(path).absolute()
    settings = pathlib.Path(settings).absolute()
    settings_dataframe = pd.read_csv(settings, sep=None, engine='python')

    # check that settings file has the correct headings
    settings_dict = settings_checker(settings_dataframe)

    # make a list of the jobs to run and the settings required and make the required folders
    list_jobs_to_run = []
    for job_entry, job_settings in settings_dict.items():
        sample_id = job_settings["sample_name"].upper()
        known_mab_name = job_settings["known_mab_name"].upper()
        run_step1 = job_settings["run_step1"]
        run_step2 = job_settings["run_step2"]
        run_step3 = job_settings["run_step3"]

        lineage = job_settings["lineage"].lower()
        time_point = job_settings["time_point"].lower()
        chain = job_settings["sonar_1_version"].lower()
        primer_name = job_settings["primer_name"].upper()

        list_jobs_to_run.append([[sample_id], [lineage, time_point, chain], [run_step1, run_step2, run_step3],
                                 [primer_name], [known_mab_name]])

        # make the folders if they don't already exist
        step_0_make_folders(path, lineage, time_point, chain, primer_name)

    # unzip_files(dir_with_files)
    unzip_files(path)

    # move files from 0_new_data, into correct directory, if necessary
    move_raw_data(path, settings_dataframe)

    command_call_processing = []
    command_call_sonar_1 = []
    command_call_sonar_2 = []

    if not list_jobs_to_run:
        print("No jobs were found in the settings file\nCheck that the file format was correct "
              "and that it contains entries")
        sys.exit("exiting")

    # counter for checking if generator is empty
    n = 0
    for job_entry in list_jobs_to_run:
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

            print("jobs to run here")
        else:
            continue

        # check that files exist in target folder
        print("Checking whether files already exist in target folder")
        search_raw_path = dir_with_raw_files.glob("*.fastq")
        search_sonar1_path = dir_with_sonar1_files.glob("*.fasta")

        # check if the generator is empty = no files in glob search
        for n, file in enumerate(search_raw_path):
            pass
        for n, file in enumerate(search_sonar1_path):
            pass

    if n == 0:
        print("No files found in target directories")
        sys.exit("exiting")
    else:
        print("files found")

    # only run sample processing if one or more files were specified
    if command_call_processing:
        step_1_run_sample_processing(command_call_processing)

    # # only run sonar1 if one or more files were specified
    # if command_call_sonar_1:
    #     step_2_run_sonar_1(command_call_sonar_1)
    #
    # # only run sonar 2 if one or more files were specified
    # if command_call_sonar_2:
    #     step_3_run_sonar_2(command_call_sonar_2)

    print("Done")


if __name__ == "__main__":
    args = docopt(__doc__, version='V0.1.0')

    main(path=args['<project_path>'], settings=args['<settings_file>'])
