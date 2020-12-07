"""Assembler wrapper for all_in.py
"""

import logging
import os
import subprocess
from typing import Dict, List, NewType, NoReturn

from tqdm.std import tqdm

FASTA = NewType('FASTA',str)
PathStr = NewType('PathStr', str)

class Assembler():
    """Base class of specific assembler
    """
    def __init__(
        self,
        fasta_name: str,
        R1_path: PathStr,
        R2_path: PathStr,
        params: List
        ) -> NoReturn:
        self.fasta_name = fasta_name
        self.R1_path = R1_path
        self.R2_path = R2_path
        self.params = params
        self.logger = logging.getLogger("all_in.assembler")
        self.logger.debug("instantiating assembler")
        self.logger.debug(
            self.__dict__
        )
    
    def assemble(self) -> FASTA:
        NotImplementedError

class Megahit(Assembler):
    """Assembler(megahit)
    """
    def __init__(self, *args, **kwargs) -> NoReturn:
        super().__init__(*args, **kwargs)
    
    def assemble(self) -> FASTA:
        # make path
        # cell_path: path to the cell
        # megahit_general_out: cell/megahit
        # megahit_out: cell/megahit/specific fasta name
        cell_path = os.path.dirname(self.R1_path)
        megahit_general_out = cell_path + '/megahit_out'
        megahit_out = megahit_general_out + '/' + self.fasta_name
        
        # Prepare directory
        if not os.path.exists(megahit_general_out):
            os.makedirs(megahit_general_out)
        
        # Construct command
        command_line = [
            "megahit",
            "-1",
            self.R1_path,
            "-2",
            self.R2_path,
            "-o",
            megahit_out,
            *self.params
        ]
        self.logger.debug("execute {}".format(command_line))
        # variables for return of proc
        msg = None
        return_code = None
        try:
            # execute
            proc = subprocess.Popen(
                command_line,
                stdin=None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )

            # get return from proc
            msg = proc.stdout.read().decode()
            return_code = proc.wait()
        except OSError as e:
            self.logger.error(msg) if msg else None
            self.logger.error("megahit returns {}".format(return_code)) if return_code else None
            self.logger.exception(e)
            res = ""
        else:
            return_code = proc.wait()
            self.logger.debug(msg)

            # load fasta
            contigs_path = megahit_out + '/final.contigs.fa'
            if os.path.exists(contigs_path):
                with open(contigs_path, 'rt') as f:
                    res = f.read()
                self.logger.debug("{} contig(s) found.".format(len(res.split('\n'))//2))
            else:
                self.logger.debug("no contig")
                res = ""
        finally:
            return res

class Skesa(Assembler):
    pass

class Spades(Assembler):
    pass

def assemble_all(
    R1_name: List[str],
    R2_name: List[str],
    destination: PathStr,
    assemble_engine: str,
    settings: Dict
    ) -> NoReturn:
    """Assemble all sequences in one cell.

    Arguments:
        R1_name(list of str): filenames in str
        R2_name(list of str): filenames in str
        destination(PathStr): path to the data folder
        assemble_engine(str): program using assembling. skesa|megahit|spades
        settings(Dict): settings load from yaml
    """
    logger = logging.getLogger("all_in.asm_function")
    logger.debug("assemble_all called.")


    # Select assembler
    try:
        selected_assembler = {
            "megahit" : Megahit,
            "skesa" : Skesa,
            "spades" : Spades
        }[assemble_engine]
    except KeyError as e:
        logger.exception(e)
        raise e
    
    # Scan "cells" directory
    cells_dir = os.path.join(destination, "cells")
    with os.scandir(cells_dir) as it:
        for entry in tqdm(list(it), desc="assemble all"):
            if entry.is_dir():
                # Construct path to the cell and each file
                path = os.path.join(cells_dir, entry.path.split('/')[-1])
                tmpR1_path = os.path.join(path, "tmp_R1.fastq")
                tmpR2_path = os.path.join(path, "tmp_R2.fastq")
                R1_fastq = [os.path.join(path, name) for name in R1_name]
                R2_fastq = [os.path.join(path, name) for name in R2_name]
                final_contigs_path = os.path.join(path, "all_contigs.fasta")

                # Make temporary fastq
                with open(tmpR1_path, 'wt') as f_tmp1, open(tmpR2_path, 'wt') as f_tmp2:
                    # read each fastq and marge into temporary fastq
                    for R1, R2 in zip(R1_fastq, R2_fastq):
                        with open(R1, 'rt') as f_r1, open(R2, 'rt') as f_r2:
                            f_tmp1.write('\n'.join(f_r1.read().split('\n')[:-1]))
                            f_tmp2.write('\n'.join(f_r2.read().split('\n')[:-1]))
                
                # construct assembler
                asm = selected_assembler(
                    fasta_name = "all",
                    R1_path = tmpR1_path,
                    R2_path = tmpR2_path,
                    params = settings[assemble_engine]
                )

                # Execute assembler
                contigs = asm.assemble()

                # Save contigs
                with open(final_contigs_path, 'wt') as f_out:
                    f_out.write(contigs)

def assemble_individually(
    R1_name: List[str],
    R2_name: List[str],
    destination: PathStr,
    assemble_engine: str,
    settings: Dict
    ) -> NoReturn:
    """Assemble sequences in one cell individually.

    Arguments:
        R1_name(list of str): filenames in str
        R2_name(list of str): filenames in str
        destination(PathStr): path to the data folder
        assemble_engine(str): program using assembling. skesa|megahit|spades
        settings(Dict): settings load from yaml
    """
    logger = logging.getLogger("all_in.asm_function")
    logger.debug("assemble_individually called.")

    # Select assembler
    try:
        selected_assembler = {
            "megahit" : Megahit,
            "skesa" : Skesa,
            "spades" : Spades
        }[assemble_engine]
    except KeyError as e:
        logger.exception(e)
        raise e

    # Scan "cells" directory
    cells_dir = os.path.join(destination, "cells")
    with os.scandir(cells_dir) as it:
        for entry in tqdm(list(it), desc="assemble individually"):
            if entry.is_dir():
                # Construct path to the cell and each file
                path = os.path.join(cells_dir, entry.path.split('/')[-1])
                R1_fastq_path = [os.path.join(path, name) for name in R1_name]
                R2_fastq_path = [os.path.join(path, name) for name in R2_name]
                common_name = [os.path.commonprefix([R1.split('.')[0], R2.split('.')[0]])[:-2] for R1,R2 in zip(R1_name, R2_name)]
                contigs_path = [os.path.join(path, name) + "_contigs.fasta" for name in common_name]
                logger.debug(
                    "path={}\nR1_fastq_path={}\nR2_fastq_path={}\ncontigs_path={}".format(
                        path, R1_fastq_path, R2_fastq_path, contigs_path
                    )
                )
                # Run assembler for each pair
                for R1, R2, name, contig in zip(R1_fastq_path, R2_fastq_path, common_name, contigs_path):
                    # Construct assembler
                    asm = selected_assembler(
                        fasta_name=name,
                        R1_path=R1,
                        R2_path=R2,
                        params=settings[assemble_engine]
                    )

                    # Execute assembler
                    contigs_fasta = asm.assemble()

                    # Save contigs
                    with open(contig, 'wt') as f_out:
                        f_out.write(contigs_fasta)