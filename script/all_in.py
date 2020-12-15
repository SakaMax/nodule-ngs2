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
import pickle
from pprint import pformat, pprint
import sys
from typing import Dict, Callable, NamedTuple, Tuple

from ruamel.yaml import YAML

import all_in_tools as tools
from all_in_tools.my_types import *

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
        required=True,
        dest="R1",
        metavar="R1.fastq"
    )
    parser.add_argument(
        "--2nd_read", "-2",
        help="The second read of paired-end sequence (fastq file)",
        action='append',
        required=True,
        dest="R2",
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
    parser.add_argument(
        "--resume_from", "-r",
        help="resume from the previous checkpoint",
        action="store",
        dest="resume",
        default=None,
        metavar="/path/to/checkpoint.json"
    )
    parser.add_argument(
        "--override_settings_by", "-o",
        help="override settings by this setting file when set with --resume_from",
        action="store",
        dest="override_yaml",
        default=None,
        metavar="/path/to/settings.yaml"
    )

    # load arguments
    # if parse_args failed, this raises SystemExit

    args = parser.parse_args()

    # return args(Namespace)
    return args

def read_settings(yaml_path : str) -> Dict:
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




class AllIn():
    """Run the script, hold status and progress of this script.

    Attributes:
        args(Namespace): arguments from command line.
        settings(dict): settings from yaml.
        logger(Logger): logger(to CLI and file)
        dir_path(dict): dict of path to directories
        fastq_path(dict): dict of path to fastq sequences
        step_counter(int): which step should this object run.
    """
    def __init__(self, args: argparse.Namespace) -> None:
        # Store arguments
        self.args = args

        # Check read
        if (args.R1 == []) or (args.R2 == []) or (len(args.R1) != len(args.R2)):
            print(
                "When the program start at the beginning, paired sequences must be given."
            )
            sys.exit(1)

        # Load settings from yaml
        self.settings = read_settings(args.settings)

        # Construct path to each dir, each fastq
        self.dir_path, self.fastq_path = construct_path(self.settings, args)

        # Initialize the logger
        self.logger = set_logger(self.settings, args)

        # Set function list
        self.step_func = [
            self._step1,
            self._step2,
            self._step3,
            self._step4,
            self._step5,
            self._step6,
            self._step7
        ]
        self.step_name = [
            "cutadapt_tag",
            "cutadapt_primer",
            "fastp",
            "demultiplex",
            "assemble_all",
            "assemble_separete",
            "blast_all"
        ]
        
        # Set counter
        self.step_counter = 0

        # Greet
        self.logger.info(
            "\n====\tThis is all_in.py for nodule NGS!\t===="
        )

        # Log current status
        self.logger.debug(pformat(self.__dict__))
    
    def workflow_generator(self) -> WorkflowFunctionInfo:
        """Yield dict about workflow function.

        returns:


        """
        start = self.step_counter
        for func, count, name in zip(
            self.step_func[start:len(self.step_func)],
            range(start, len(self.step_func)),
            self.step_name[start:len(self.step_func)]
        ):
            print("before yield {}".format(name))
            # yield {
            #     "function": func,
            #     "name" : name,
            #     "step": count+1
            # }
            yield WorkflowFunctionInfo(
                function=func,
                name=name,
                step=count+1,
            )
            print("after yield {}".format(name))
            self.step_counter = count + 1
            # try:
            #     self.logger.info("====STEP {}: {}====".format(count+1,name))
            #     func()
            # except Exception as e:
            #     self.logger.exception(e)
            #     raise e
            # else:
            #     self.logger.info("===={} end.====".format(name))
            #     self.step_counter = count + 1
            #     yield self.step_counter

    def _step1(self) -> None:
        """Cutadapt for tag recognition
        
        """
        tools.cutadapt.cutadapt_tag(
            R1_fastq=self.fastq_path["raw"]["R1"],
            R2_fastq=self.fastq_path["raw"]["R2"],
            forward_tag=self.dir_path["tag"][0],
            reverse_tag=self.dir_path["tag"][1],
            destination=self.dir_path["destination"]
        )

    def _step2(self) -> None:
        """Cutadapt for common primer
        """
        tools.cutadapt.cutadapt_primer(
            R1_fastq=self.fastq_path["tag_removed"]["R1"],
            R2_fastq=self.fastq_path["tag_removed"]["R2"],
            forward_primer=self.dir_path["primer"][0],
            reverse_primer=self.dir_path["primer"][1],
            destination=self.dir_path["destination"]
        )

    def _step3(self) -> None:
        """Quality filtering by fastp
        """
        tools.fastp.fastp(
            R1_fastq=self.fastq_path["primer_removed"]["R1"],
            R2_fastq=self.fastq_path["primer_removed"]["R2"],
            destination=self.dir_path["destination"],
            report_dest=self.dir_path["report_dest"],
            settings=self.settings
        )

    def _step4(self) -> None:
        """Demultiplex
        """
        tools.demultiplex.demultiplex(
            R1_fastq=self.fastq_path["fastp"]["R1"],
            R2_fastq=self.fastq_path["fastp"]["R2"],
            cells_json=self.settings["data"]["cells_json"],
            destination=self.dir_path["destination"]
        )

    def _step5(self) -> None:
        """Assemble (Using all sequences at one time)
        """
        tools.assemble.assemble_all(
            R1_name = [r1.split('/')[-1] for r1 in self.fastq_path["fastp"]["R1"]],
            R2_name = [r2.split('/')[-1] for r2 in self.fastq_path["fastp"]["R2"]],
            destination = self.dir_path["destination"],
            assemble_engine = self.args.engine,
            settings=self.settings
        )

    def _step6(self) -> None:
        tools.assemble.assemble_individually(
            R1_name = [r1.split('/')[-1] for r1 in self.fastq_path["fastp"]["R1"]],
            R2_name = [r2.split('/')[-1] for r2 in self.fastq_path["fastp"]["R2"]],
            destination = self.dir_path["destination"],
            assemble_engine = self.args.engine,
            settings=self.settings
        )
    
    def _step7(self) -> None:
        """Blast search(all)
        """
        self.blast_all: dict[str: BlastResultInfo] = tools.blast.blast_all(
            destination = self.dir_path["destination"],
            settings=self.settings
        )

class AllInManager():
    """Load/Save AllIn and run AllIn's workflow

    Attributes:
        args(Namespace): arguments from stdin
        _all_in(AllIn): actual scripts and its status
    """
    def __init__(self, args: dict) -> None:
        self.args = args

        if args.resume is not None:
            # resume point
            with open(args.resume, 'rb') as f:
                self._all_in = pickle.load(f)

            # If override flag is on, reload settings
            if args.override_yaml is not None:
                self._all_in.settings = read_settings(args.override_yaml)

            # re-setup logger
            self._all_in.logger = set_logger(
                self._all_in.settings, self._all_in.args
            )
            self._all_in.logger.info("All_In resumed")
        else:
            # new start
            self._all_in = AllIn(args)
            self._all_in.logger.info("All_In intialized")

    def run(self) -> None:
        """Run the script
        """
        # step_info is WorkflowFunctionInfo.
        # for each step,
        # - save current status
        # - log step No. and name
        # - run the step
        # - report exception if it occurs
        # - log end of the step
        for step_info in self._all_in.workflow_generator():
            self._create_checkpoint("before_" + step_info.name)
            try:
                self._all_in.logger.info(
                    "====step {}: {}====".format(
                        step_info.step, step_info.name)
                    )
                step_info.function()
            except Exception as e:
                self._all_in.logger.exception(e)
            else:
                self._all_in.logger.info("step {} end".format(step_info.step))
        else:
            # Create checkpoint for test
            self._create_checkpoint("after_all")

    def _create_checkpoint(self, name) -> None:
        """Save self._all_in to 'destination/name.checkpoint'

        Arguments:
            name(str): name of the checkpoint
        """
        if not os.path.exists(self._all_in.dir_path["destination"]):
            os.makedirs(self._all_in.dir_path["destination"])
        
        chk_path = os.path.join(
            self._all_in.dir_path["destination"], name + ".checkpoint"
        )
        with open(chk_path, 'wb') as f:
            pickle.dump(self._all_in, f)

if __name__ == "__main__":
    # read args
    args = get_args()

    # Get AllInManager:
    # Inside this, manager creates AllIn (if there is no checkpoint in arguments) or
    # loads AllIn (otherwise)
    manager = AllInManager(args)

    # Start actual script
    # Before each step, the manager automatically create checkpoints
    manager.run()