import logging
import os
from os import PathLike
import subprocess
from typing import List, NoReturn

def cutadapt_tag(
    R1_fastq: List[PathLike],
    R2_fastq: List[PathLike],
    forward_tag: PathLike,
    reverse_tag: PathLike,
    destination: PathLike,
) -> NoReturn:
    """Run cutadapt to recognize tag

    This function execute cutadapt and store tag-removed fastq in destination/tag_removed.

    Arguments:
        R1_fastq(list of PathLike): path to the raw R1 sequence
        R2_fastq(list of PathLike): path to the raw R2 sequence
        forward_tag(PathLike): path to the forward(R1) tag 
        reverse_tag(PathLike): path to the reverse(R2) tag
        destination(PathLike): path to the data folder
    """
    logger = logging.getLogger("all_in.cutadapt")
    logger.debug("cutadapt_tag called.")

    # Prepare destination directory
    tag_removed_path = os.path.join(
        destination, "tag_removed"
    )
    if not os.path.exists(tag_removed_path):
        os.makedirs(tag_removed_path)

    # run cutadapt
    for R1, R2 in zip(R1_fastq, R2_fastq):
        # get filename
        R1_name = R1.split('/')[-1].split('.')[0]
        R2_name = R2.split('/')[-1].split('.')[0]
        # Create command for cutadapt
        command_line = [
            "python3",
            "-m",
            "cutadapt",
            "--no-indels",
            "--discard-untrimmed",
            "-g",
            "file:{}".format(forward_tag),
            "-G",
            "file:{}".format(reverse_tag),
            "-y",
            r" {name}",
            "-o",
            "{}/tag_removed/{}.fastq".format(destination, R1_name),
            "-p",
            "{}/tag_removed/{}.fastq".format(destination, R2_name),
            R1,
            R2
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
            logger.error("cutadapt returns {}".format(return_code)) if return_code else None
            logger.exception(e)
            raise e
        else:
            return_code = proc.wait()
            logger.debug(msg)

def cutadapt_primer(
    R1_fastq: List[PathLike],
    R2_fastq: List[PathLike],
    forward_primer: PathLike,
    reverse_primer: PathLike,
    destination: PathLike,
) -> NoReturn:
    """Run cutadapt to remove common primers

    This function execute cutadapt and store primer-removed fastq in destination/primer_removed.

    Arguments:
        R1_fastq(PathLike): path to the tag_removed R1 sequence
        R2_fastq(PathLike): path to the tag_removed R2 sequence
        forward_primer(PathLike): path to the forward(R1) tag 
        reverse_primer(PathLike): path to the reverse(R2) tag
        destination(PathLike): path to the data folder
    """
    logger = logging.getLogger("all_in.cutadapt")
    logger.debug("cutadapt_primer called.")

    # Prepare destination directory
    primer_removed_path = os.path.join(
        destination, "primer_removed"
    )
    if not os.path.exists(primer_removed_path):
        os.makedirs(primer_removed_path)

    # run cutadapt
    for R1, R2 in zip(R1_fastq, R2_fastq):
        # get filename
        R1_name = R1.split('/')[-1].split('.')[0]
        R2_name = R2.split('/')[-1].split('.')[0]
        # Create command for cutadapt
        command_line = [
            "python3",
            "-m",
            "cutadapt",
            "--no-indels",
            "--discard-untrimmed",
            "-g",
            "file:{}".format(forward_primer),
            "-G",
            "file:{}".format(reverse_primer),
            "-o",
            "{}/primer_removed/{}.fastq".format(destination, R1_name),
            "-p",
            "{}/primer_removed/{}.fastq".format(destination, R2_name),
            R1,
            R2
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
            logger.error("cutadapt returns {}".format(return_code)) if return_code else None
            logger.exception(e)
            raise e
        else:
            return_code = proc.wait()
            logger.debug(msg)