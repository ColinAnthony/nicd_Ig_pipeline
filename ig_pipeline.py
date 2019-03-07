#!/anaconda/bin/python3
"""This is a wrapper to automate the NGS nAb processing and analysis pipeline.

    Usage:
         ig_pipeline.py (-p <project_path>) (-s <settings_file>)
         ig_pipeline.py (-p <project_path>) (-s <settings_file>) [-r <rerun_step>]
         ig_pipeline.py -h
         ig_pipeline.py -v

    Options:
         -p --path         The path to the project folder, where the folders will be created
         -s --settings     The path and file name of the settings csv file
         -r --run_step
         -v --version      Show the script version number
         -h --help         Show this screen.
"""
import sys
from docopt import docopt
import pandas as pd
import pathlib
import shutil
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
        try:
            parts = sample_id.split("_")
            if len(parts) != 5:
                print(f"your sample name was {sample_id}\n")
                sys.exit()
            elif parts[0][:3].upper() == "CAP" and len(parts[0]) != 6:
                print(f"your sample name was {sample_id}\n The PID must be zero padded: eg: CAP008, not CAP8" )
                sys.exit()
            elif parts[0][:3].upper() == "CAP" and len(parts[1]) != 4:
                print(f"your sample name was {sample_id}\n The visit code must be in the 2000 format, not P2V0" )
                sys.exit()
            elif parts[0][:3].upper() == "CAP" and len(parts[2]) != 6 and parts[2][-3:].lower() != "wpi":
                print(f"your sample name was {sample_id}\n The wpi code must be zero padded and end in wpi eg: 080wpi" )
                sys.exit()

        except:
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


def move_raw_data(path, settings_dataframe):
    """
    function to move raw data to the correct folder within the project
    :param path: path to the project folder
    :param settings_dataframe: the dataframe object form the settings file
    :return: a list of all the files that were moved
    """
    raw_files = pathlib.Path(path, "0_new_data").glob("*R1*.fastq*")
    all_destinations = []
    if raw_files is not None:
        for fileR1 in raw_files:
            name = fileR1.stem

            name = "_".join(name.split("_")[:5]).strip()
            print(name)

            file_extensions = fileR1.suffixes
            file_ending = ''
            for end in file_extensions:
                file_ending = file_ending + end

            fields = settings_dataframe.loc[settings_dataframe['sample_name'] == name].head(1)

            lineage = fields.iloc[0]["lineage"].lower()
            chain = fields.iloc[0]["sonar_1_version"].lower()
            time_point = fields.iloc[0]["time_point"].lower()
            destination_R1 = pathlib.Path(path, lineage, time_point, chain, name + "_R1_" + file_ending)
            print(destination_R1)
            all_destinations.append(destination_R1)
            if pathlib.Path(destination_R1).is_file():
                print("already exists")
            else:
                shutil.move(str(fileR1), str(destination_R1))
                fileR2 = str(fileR1).replace("_R1_", "_R2_")
                destination_R2 = str(destination_R1).replace("_R1_", "_R2_")
                all_destinations.append(destination_R2)
                shutil.move(str(fileR2), str(destination_R2))

    return sorted(all_destinations)


def unzip_files(all_destinations):
    for file in all_destinations:
        suffix = file.suffixe
        if suffix == "zip":
            cmd = f"unzip {file}"
            subprocess.call(cmd, shell=True)
        elif suffix == "gz":
            cmd = f"gunzip {file}"
            subprocess.call(cmd, shell=True)
        elif suffix == "fastq":
            print("the file is already unziped")


def step_1_run_sample_processing():

    pass


def step_2_run_sonar_1():
    pass


def step_3_run_sonar_2():
    pass


def main(path, settings):
    """
    create folder structure for a sequencing project
    :param your_path: (str) path to where the folders should be made
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

    for job_entry, job_settings in settings_dict.items():
        sample_id = job_settings["sample_name"].upper()
        lineage = job_settings["lineage"].lower()
        time_point = job_settings["time_point"].lower()
        chain = job_settings["sonar_1_version"].lower()
        primer_name = job_settings["primer_name"].upper()
        known_mab_name = job_settings["known_mab_name"].upper()

        run_step1 = job_settings["run_step1"]
        run_step2 = job_settings["run_step2"]
        run_step3 = job_settings["run_step3"]

        # make the folders if they don't already exist
        step_0_make_folders(path, lineage, time_point, chain, primer_name)

        # move files from 0_new_data, into correct directory, if necessary
        all_destinations = move_raw_data(path, settings_dataframe)
        print(all_destinations)
        # unzip_files(all_destinations)

        step_1_run_sample_processing()


        step_2_run_sonar_1()



    print("Done")


if __name__ == "__main__":
    args = docopt(__doc__, version='V0.1.0')
    # print(args)

    main(path=args['<project_path>'], settings=args['<settings_file>'])

    # main(docopt(arguments))

    # print(args)
    # main(path, settings)

    # parser = argparse.ArgumentParser(description='This is a wrapper to automate the NGS nAb processing and analysis '
    #                                              'pipeline')
    # parser.add_argument('-p', '--path', type=str,
    #                     help='The path to the project folder, where the folders will be created', required=True)
    # parser.add_argument('-s', '--settings', type=str,
    #                     help='the path and file name of the settings csv file', required=True)
    #
    # args = parser.parse_args()
    # path = args.path
    # settings = args.settings
    #
    # main(path, settings)
