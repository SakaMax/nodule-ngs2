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
import sys
from typing import NoReturn

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
    try:
        with open(yaml_path, 'rt') as f:
            user_settings = yaml.load(f)
    except FileNotFoundError:
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

    # construct path:
    # path to repository, tag, primer, cell list, data destination
    repo_dir = '/'.join(
        os.path.dirname(
            os.path.abspath(__file__)
        ).split('/')[:-1]
    )
    tag_path = [
        os.path.join(repo_dir, tag) for tag in settings['data']['tag']
    ]
    primer_path = [
        os.path.join(repo_dir, primer) for primer in settings['data']['primer']
    ]
    cells_path = os.path.join(repo_dir, settings['data']['cells_json'])
    destination_path = os.path.join(
        repo_dir,
        settings['data']['destination'],
        datetime.now().strftime(settings['data']['datetime_format'])
    )
    logger.debug(
        """
        repo        ->{}
        tag         ->{}
        primer      ->{}
        cells       ->{}
        destination ->{}""".format(
            repo_dir, tag_path, primer_path, cells_path, destination_path
        )
    )

    # STEP 1
    # Recognition and removal of tag
    logger.info("STEP 1: cutadapt (for tag)")
    tools.cutadapt.cutadapt_tag(
        args.R1, args.R2, tag_path[0], tag_path[1],
        destination_path, settings
    )