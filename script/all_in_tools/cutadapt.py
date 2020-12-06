import logging
import os
from os import PathLike
import subprocess
from typing import Dict, NoReturn

def cutadapt_tag(
    R1_fastq: PathLike,
    R2_fastq: PathLike,
    forward_tag: PathLike,
    reverse_tag: PathLike,
    destination: PathLike,
) -> NoReturn:
    """Run cutadapt to recognize tag

    This function execute cutadapt and store tag-removed fastq in destination/tag_removed.

    Arguments:
        R1_fastq(PathLike): path to the raw R1 sequence
        R2_fastq(PathLike): path to the raw R2 sequence
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


        try:
            proc = subprocess.Popen(
            command_line,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
            )
            msg = proc.stdout.read().decode()
        except OSError as e:
            return_code = proc.wait()
            logger.error(msg)
            logger.error("cutadapt returns {}".format(return_code))
            logger.exception(e)
            raise e
        else:
            return_code = proc.wait()
            logger.info(msg)

def cutadapt_primer(
    R1_fastq: PathLike,
    R2_fastq: PathLike,
    forward_primer: PathLike,
    reverse_primer: PathLike,
    destination: PathLike,
) -> NoReturn:
    """Run cutadapt to remove common primers

    This function execute cutadapt and store primer-removed fastq in destination/primer_removed.

    Arguments:
        R1_fastq(PathLike): path to the raw R1 sequence
        R2_fastq(PathLike): path to the raw R2 sequence
        forward_tag(PathLike): path to the forward(R1) tag 
        reverse_tag(PathLike): path to the reverse(R2) tag
        destination(PathLike): path to the data folder
    """