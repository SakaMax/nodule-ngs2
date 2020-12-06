"""Demultiplex paired-end fastq file after fastp.

Note:
    dependencies(numpy, pandas)
    This script can't handle the interleave format.
"""
from datetime import datetime
import json
import logging
from os import PathLike
import os
import sys
from typing import List, NoReturn

import numpy as np
import pandas as pd

def demultiplex(
    R1_fastq: List[PathLike],
    R2_fastq: List[PathLike],
    cells_json: PathLike,
    destination: PathLike,
) -> NoReturn:
    """Demultiplex fastq file (after fastp)

    Arguments:
        R1_fastq(PathLike): path to the after-fastp R1 sequence
        R2_fastq(PathLike): path to the after-fastp R2 sequence
        destination(PathLike): path to the data folder
        cells_json(PathLike): path to ``cells.json``
    """

    logger = logging.getLogger("all_in.demultiplex")
    logger.debug("demultiplex called")

    # Load cell names
    try:
        with open(cells_json, 'rt') as f:
            cells = json.load(f)
    except FileNotFoundError as e:
        logger.exception(e)
        logger.fatal("Can't open cell list. Abort.")
        sys.exit(1)
    
    for R1, R2 in zip(R1_fastq, R2_fastq):
        # Read two fastq file
        # Fastq files generated by fastp has one blank line at the end,
        # so we need all lines in the fastq file except the last line.
        # Because of this, split('\n')[:-1] is used here.
        with open(R1, 'rt') as f_r1, open(R2, 'rt') as f_r2:
            fastq_r1 = f_r1.read().split('\n')[:-1]
            fastq_r2 = f_r2.read().split('\n')[:-1]

        # ensure both R1 and R2 have same length
        try:
            assert len(fastq_r1) == len(fastq_r2)
        except AssertionError as e:
            logger.fatal(
                "R1({} lines) != R2({} lines)".format(len(fastq_r1), len(fastq_r2))
                )
            logger.exception(e)
            exit(1)
        else:
            fastq_length = len(fastq_r1)
            logger.debug(
                "fastq length assertion OK.(length={})".format(fastq_length)
                )
        
        # Convert fastq into DataFrame via dict
        dict_r1 = {
            "R1_fullname" : [fastq_r1[i] for i in range(0,fastq_length,4)],
            "R1_name"     : [fastq_r1[i].split(" ")[-1] for i in range(0,fastq_length,4)],
            "R1_sequence" : [fastq_r1[i+1] for i in range(0,fastq_length,4)],
            "R1_quality"  : [fastq_r1[i+3] for i in range(0,fastq_length,4)]
        }
        dict_r2 = {
            "R2_fullname" : [fastq_r2[i] for i in range(0,fastq_length,4)],
            "R2_name"     : [fastq_r2[i].split(" ")[-1] for i in range(0,fastq_length,4)],
            "R2_sequence" : [fastq_r2[i+1] for i in range(0,fastq_length,4)],
            "R2_quality"  : [fastq_r2[i+3] for i in range(0,fastq_length,4)]
        }
        sequence_df = pd.concat(
            [pd.DataFrame(dict_r1), pd.DataFrame(dict_r2)],
            axis=1
        )
        #sequence_df.to_csv("test_df.csv")
        logger.debug(sequence_df.head())

        # split sequence_df into cell dfs
        # {
        #   "1A01" : dataframe contains sequences from 1A01,
        #   "1A02" : dataframe contains sequences from 1A02,
        #   ...
        # }
        cell_dfs = {
            k: sequence_df.query(
                " or".join(
                    ["(R1_name == '{}' and R2_name == '{}')".format(pair[0],pair[1]) \
                        for pair in cells[k]]
                    )
            ) for k in cells.keys()
        }

        # Count sequences in each cell
        plate = [
            pd.DataFrame(
                np.zeros((8,12),dtype=int), 
                index=["A","B", "C", "D", "E", "F", "G", "H"], 
                columns=["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
            ) for _ in range(0,6)
        ]
        empty_cells = list()
        for cell, df in cell_dfs.items():
            plate_no = int(cell[0])-1
            row = cell[1]
            col = cell[2:4]
            plate[plate_no].at[row, col] = len(df)
            if len(df) == 0:
                empty_cells.append(cell)

        # Print plate shape
        for i in range(0,6):
            print("\t==== Reads in plate No. {} ====".format(i+1))
            print(plate[i])
        else:
            print("Empty cells : {}".format(empty_cells))

        logger.info("empty cells : {} out of {}".format(len(empty_cells), len(cell_dfs)))
        # Save data
        # destination/
        #   |-cells/
        #       |-1A01/
        #       |   |-R1.fastq
        #       |   |-R2.fastq
        #       |-1A02/
        #       ...
        for cell, df in cell_dfs.items():
            path = destination + '/cells/' + cell
            os.makedirs(path) if not os.path.exists(path) else None

            with open(path + '/R1.fastq', 'wt') as f_r1, \
                open(path + '/R2.fastq', 'wt') as f_r2:

                for _, row in df.iterrows():
                    # Reproduce the sequence in fastq format
                    r1_seq = '\n'.join(
                        [row.R1_fullname, row.R1_sequence, '+', row.R1_quality]
                    )
                    r2_seq = '\n'.join(
                        [row.R2_fullname, row.R2_sequence, '+', row.R2_quality]
                    )
                    
                    # Write sequence
                    f_r1.write(r1_seq + '\n')
                    f_r2.write(r2_seq + '\n')

