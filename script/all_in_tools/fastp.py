import logging
import os
import shutil
import subprocess
from typing import Dict, List , NewType

PathStr = NewType('PathStr', str)

def fastp(
    R1_fastq: List[PathStr],
    R2_fastq: List[PathStr],
    destination: PathStr,
    report_dest: PathStr,
    settings: Dict,
) -> None:
    """Run fastp to quality filtering

    Arguments:
        R1_fastq(PathStr): path to the primer_removed R1 sequence
        R2_fastq(PathStr): path to the primer_removed R2 sequence
        destination(PathStr): path to the data folder
        report_dest(PathStr): path to the nginx html folder
        settings(dict): settings from yaml
    """

    logger = logging.getLogger("all_in.fastp")
    logger.debug("fastp called")

    # Prepare destination directory
    fastp_path = os.path.join(
        destination, "fastp"
    )
    if not os.path.exists(fastp_path):
        os.makedirs(fastp_path)
    if not os.path.exists(report_dest):
        os.makedirs(report_dest)

    # Run fastp
    for R1, R2 in zip(R1_fastq, R2_fastq):
        # Get filename
        R1_name = R1.split('/')[-1].split('.')[0]
        R2_name = R2.split('/')[-1].split('.')[0]
        # Get common name
        common_name = os.path.commonprefix([R1_name, R2_name])

        # Create command for fastp
        command_line = [
            "fastp",
            "-i",
            R1,
            "-I",
            R2,
            "-o",
            os.path.join(fastp_path, R1_name) + ".fastq",
            "-O",
            os.path.join(fastp_path, R2_name) + ".fastq",
            "-h",
            os.path.join(fastp_path, common_name) + ".html",
            *settings["fastp"]
        ]
        logger.debug("execute {}".format(command_line))

        # variables for return of proc
        msg = None
        return_code = None

        try:
            # execute proc
            proc = subprocess.Popen(
            command_line,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
            )

            # get return from proc
            msg = proc.stdout.read().decode()
            return_code = proc.wait()
        except OSError as e:
            logger.error(msg) if msg else None
            logger.error("fastp returns {}".format(return_code)) if return_code else None
            logger.exception(e)
            raise e
        else:
            return_code = proc.wait()
            logger.info(msg)
            shutil.copy(
                os.path.join(fastp_path, common_name) + ".html",
                report_dest
            )