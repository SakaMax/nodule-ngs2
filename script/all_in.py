"""Entry point of workflow: replacement for all_in_one.sh

Features...
    * Tag & Adapter trimming
    * Quality filtering
    * De novo assembling
    * Homology searching

This script reads settings from the  ``settings.yaml`` .
And, of course, fastq files are needed. 

Example:
    If you have a paired-read fastq (R1.fastq and R2.fastq) sequence and use SKESA,
    type (in shell)::

        $ python all_in.py -1 R1.fastq -2 R2.fastq -a skesa

    When you assemble using megahit and apply some parameters to megahit::

        $ python all_in.py -1 R1.fastq -2 R2.fastq -a megahit -s ./settings.yaml

    For more infomation::

        $ python all_in.py --help
    
    or check ``settings.yaml`` .

    
"""

import argparse
from datetime import datetime
import logging
import os
from pprint import pformat
import sys
from typing import NoReturn, Tuple

from ruamel.yaml import YAML

import all_in_tools as tools

def get_args() -> argparse.Namespace:
    """Read given arguments.

    returns:
        argparse.Namespace : given arguments from commandline.
    """

    # Create new ArgumentParser
    parser = argparse.ArgumentParser()

    # Add arguments
    parser.add_argument(
        "--1st-read","-1",
        help="The first read of paired-end sequence (fastq file)",
        action='append',
        dest="R1",
        required=True,
        metavar="R1.fastq"
    )
    parser.add_argument(
        "--2nd_read", "-2",
        help="The second read of paired-end sequence (fastq file)",
        action='append',
        dest="R2",
        required=True,
        metavar="R2.fastq"
    )
    parser.add_argument(
        "--assembler", "-a",
        help="assemble method",
        action="store",
        dest="engine",
        choices=["skesa", "megahit", "spades"],
        default="megahit",
    )
    parser.add_argument(
        "--settings", "-s",
        help="where the settings.yaml is",
        action="store",
        default=None,
        metavar="/path/to/settings.yaml"
    )

    # load arguments
    # if parse_args failed, this raises SystemExit

    args = parser.parse_args()

    # return args(Namespace)
    return args

def read_settings(yaml_path : str) -> dict:
    """read ``settings.yaml`` from yaml_path.

    arguments:
        yaml_path(str) : where the yaml file is.

    return:
        dict: yaml data.
    """
    # construct yaml parser
    yaml = YAML(typ='safe')

    # load yaml file.
    if yaml_path is not None:
        try:
            with open(yaml_path, 'rt') as f:
                user_settings = yaml.load(f)
        except FileNotFoundError:
            user_settings = dict()
    else:
        user_settings = dict()

    # load default yaml
    default_yaml_path = os.path.join(os.path.dirname(__file__),'default.yaml')
    try:
        with open(default_yaml_path, 'rt') as f:
            default_settings = yaml.load(f)
    except FileNotFoundError:
        print("{} not found. abort.".format(default_yaml_path))
        sys.exit(1)

    # marge user settings and default settings
    # if user settings exist, override default settings
    settings = {
        key: user_settings[key] if key in user_settings else default_settings[key] \
            for key in default_settings.keys()
    }

    return settings

def set_logger(settings: dict, args: argparse.Namespace) -> logging.Logger:
    """set logging settings (logging.basicConfig)

    arguments:
        settings(dict): settings from yaml.
        args(argparse.Namespace): arguments.
    returns:
        None
    """
    # Log levels
    levels = {
        "DEBUG" : logging.DEBUG,
        "INFO" : logging.INFO,
        "WARNING" : logging.WARNING,
        "ERROR" : logging.ERROR
    }

    # Prepare log directory
    if not os.path.exists(settings["logging_path"]):
        os.makedirs(
            os.path.dirname(settings["logging_path"])
        )

    # get root logger
    logger = logging.getLogger("all_in")
    logger.setLevel(logging.DEBUG)

    # create log handler
    # fh : filehandler
    # ch : consolehandler
    fh = logging.FileHandler(
        os.path.join(settings["logging_path"], settings["filename"]),
        mode=settings["logging_config"]["filemode"]
    )
    ch = logging.StreamHandler()

    # set logging level
    fh.setLevel(levels[settings["logfile_level"]])
    ch.setLevel(levels[settings["console_level"]])
    
    # set formatter
    formatter = logging.Formatter(
        fmt=settings["logging_config"]["format"],
        datefmt=settings["logging_config"]["datefmt"]
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add handlers to logger
    logger.addHandler(fh)
    logger.addHandler(ch)


    logger.debug("logging configured.")
    logger.debug("arguments : \n\t{}".format(args))
    logger.debug("settings : \n\t{}".format(settings))

    return logger

def construct_path(settings: dict, args: argparse.Namespace) -> Tuple[dict, dict]:
    """Construct path to directories, fastq files.

    arguments:
        settings(dict): settings from yaml.
        args(argparse.Namespace): arguments.
    returns:
        tuple(dict,dict): (dict of directory path, dict of fastq path)  
    """

    repo_dir = '/'.join(
        os.path.dirname(
            os.path.abspath(__file__)
        ).split('/')[:-1]
    )
    time_prefix = datetime.now().strftime(settings['data']['datetime_format'])

    # Path to directory
    dir_dict = {
        "repo" : repo_dir,
        "tag" : [
        os.path.join(repo_dir, tag) for tag in settings['data']['tag']
        ],
        "primer" : [
        os.path.join(repo_dir, primer) for primer in settings['data']['primer']
        ],
        "cells" : os.path.join(repo_dir, settings['data']['cells_json']),
        "destination" : os.path.join(
            repo_dir,
            settings['data']['destination'],
            time_prefix
        ),
        "report_dest" : os.path.join(
            settings["fastp_output"],
            time_prefix
        )
    }

    # Path to fastq
    fastq_dict = {
        "raw" : {
            "R1" : args.R1, "R2" : args.R2
        },
        "tag_removed" : {
            "R1" : [
                os.path.join(
                    dir_dict["destination"], "tag_removed", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R1
            ],
            "R2" : [
                os.path.join(
                    dir_dict["destination"], "tag_removed", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R2
            ]
        },
        "primer_removed" : {
            "R1" : [
                os.path.join(
                    dir_dict["destination"], "primer_removed", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R1
            ],
            "R2" : [
                os.path.join(
                    dir_dict["destination"], "primer_removed", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R2
            ]
        },
        "fastp" : {
            "R1" : [
                os.path.join(
                    dir_dict["destination"], "fastp", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R1
            ],
            "R2" : [
                os.path.join(
                    dir_dict["destination"], "fastp", '.'.join(f.split('/')[-1].split('.')[0:2])
                    ) for f in args.R2
            ]
        }
    }

    return dir_dict, fastq_dict

if __name__ == "__main__":
    # read args
    args = get_args()

    # read settings
    settings = read_settings(args.settings)

    # set up logging
    logger = set_logger(settings, args)

    # greet
    logger.info(
        "\n====\tThis is all_in.py for nodule NGS!\t===="
    )

    # construct path
    dir_path, fastq_path = construct_path(settings, args)

    logger.debug(pformat(dir_path))
    logger.debug(pformat(fastq_path))

    # STEP 1
    # Recognition and removal of tag
    logger.info("STEP 1: cutadapt (for tag)")
    tools.cutadapt.cutadapt_tag(
        R1_fastq=fastq_path["raw"]["R1"],
        R2_fastq=fastq_path["raw"]["R2"],
        forward_tag=dir_path["tag"][0],
        reverse_tag=dir_path["tag"][1],
        destination=dir_path["destination"]
    )

    # STEP 2
    # Recognition and removal of primer
    logger.info("STEP 2: cutadapt (for primer)")
    tools.cutadapt.cutadapt_primer(
        R1_fastq=fastq_path["tag_removed"]["R1"],
        R2_fastq=fastq_path["tag_removed"]["R2"],
        forward_primer=dir_path["primer"][0],
        reverse_primer=dir_path["primer"][1],
        destination=dir_path["destination"]
    )

    # STEP 3
    # Quality filtering (fastp)
    logger.info("STEP 3: fastp")
    tools.fastp.fastp(
        R1_fastq=fastq_path["primer_removed"]["R1"],
        R2_fastq=fastq_path["primer_removed"]["R2"],
        destination=dir_path["destination"],
        report_dest=dir_path["report_dest"],
        settings=settings
    )

    # STEP 4
    # Demultiplex
    logger.info("STEP 4: Demultiplex")
    tools.demultiplex.demultiplex(
        R1_fastq=fastq_path["fastp"]["R1"],
        R2_fastq=fastq_path["fastp"]["R2"],
        cells_json=settings["data"]["cells_json"],
        destination=dir_path["destination"]
    )

    # STEP 5
    # Assemble
    logger.info("STEP 5: Assemble")
    tools.assemble.assemble_all(
        R1_name = [r1.split('/')[-1] for r1 in fastq_path["fastp"]["R1"]],
        R2_name = [r2.split('/')[-1] for r2 in fastq_path["fastp"]["R2"]],
        destination = dir_path["destination"],
        assemble_engine = args.engine,
        settings=settings
    )
    tools.assemble.assemble_individually(
        R1_name = [r1.split('/')[-1] for r1 in fastq_path["fastp"]["R1"]],
        R2_name = [r2.split('/')[-1] for r2 in fastq_path["fastp"]["R2"]],
        destination = dir_path["destination"],
        assemble_engine = args.engine,
        settings=settings
    )