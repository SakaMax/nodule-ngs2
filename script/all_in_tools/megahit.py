""" SKESA wrapper for NGS workflow

Example:
    To run this, give one argument: path/to/cells_folder
    When you assemble set of sequences in ./data/2020_1201_1020/cells, type
    ::
        $ python skesa.py ./data/2020_1201_1020/cells

    Then the script automatically find folders like 1A01, 1A02, ...
"""

import argparse
from datetime import datetime
import logging
from logging import exception, log
import os
import subprocess
from subprocess import PIPE
import sys

from tqdm import tqdm

def get_args() -> argparse.Namespace:
    """Read given arguments.

    returns:
        argparse.Namespace : given arguments from commandline.
    """

    # create new ArgumentParser
    parser = argparse.ArgumentParser()

    # add 3 arguments:
    # R1 fastq file, R2 fastq file,
    # destination folder
    # all arguments are positional.
    parser.add_argument("cells_path", help="path to the 'cells' folder created after running demultiplex.py.")

    # load arguments
    # if parse_args failed, this raises SystemExit
    try:
        args = parser.parse_args()
    except SystemExit as sys_exit:
        logging.error("Invalid argv{}".format(sys.argv))
        logging.exception(sys_exit)
        exit(sys_exit.code)
    else:
        logging.info(args)
    
    # return args(Namespace)
    return args

if __name__ == "__main__":
        # Activate logging
    logging.basicConfig(
        filename="megahit.log",
        format="%(levelname)s:%(message)s",
        level=logging.DEBUG
    )

    try:
        # Record time
        logging.info("---- start megahit.py {} ----".format(datetime.now()))

        # Get arguments
        args = get_args()

        # Scan 'cells' directory
        with os.scandir(args.cells_path) as it:
            for entry in tqdm(list(it)):
                if entry.is_dir():
                    # Construct the path to the cell
                    path = args.cells_path + '/' + entry.path.split('/')[-1]
                    # Construct command
                    command_line = [
                        "megahit",
                        "--keep-tmp-files",
                        "-1",
                        "{}/R1.fastq".format(path),
                        "-2",
                        "{}/R2.fastq".format(path),
                        "-o",
                        "{}".format(path + "/megahit_out")
                    ]
                    logging.debug("execute {}".format(command_line))
                    # Run SKESA with logfile
                    with open(path + "/megahit.log", 'wt') as f:
                        try:
                            # Call new process
                            proc = subprocess.Popen(
                                command_line,
                                stdin=None,
                                stdout=f,
                                stderr=f
                            )
                        except OSError as e:
                            logging.exception(e)
                            continue
                        else:
                            return_code = proc.wait()
                            if return_code != 0:
                                # something wrong
                                logging.warning(
                                    "{} returns {}".format(command_line, return_code)
                                    )
                            else:
                                logging.debug(
                                    "{} ... OK.".format(entry.path.split('/')[-1])
                                )
                    
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
    else:
        logging.info("Everything seems OK.")
    finally:
        logging.info("---- end {} ----\n".format(datetime.now()))