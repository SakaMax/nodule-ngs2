"""Blastn wrapper for all_in.py
"""

import logging
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Dict, List, NamedTuple, NewType, Optional, Union
import sys
from numpy.lib.arraysetops import isin

import pandas as pd
from tqdm.std import tqdm

# Type alias
FASTA = str
PathStr = str
BlastResultDF = pd.DataFrame

class BlastResultInfo(NamedTuple):
    """Blast result and additional infomation.

    Attributes:
        result(BlastResultDF): search result.
        query_file(PathStr): query fasta file.
        query_name(str): query name i.e. >THIS
        query_seq(str): DNA sequence used.
    """
    result: BlastResultDF
    query_file: PathStr
    query_name: str
    query_seq: str

class Blastn():
    """Blastn search engine.

    Attributes:
        settings(dict): settings from yaml.
        params(list): parameters passed to blastn.

    Note:
        To execute search, use "blast_search."

        BlastResult can be None if there is no match.
        (This is why the return of search functions is "Optional.")
    """
    def __init__(
        self,
        settings: Dict,
        params: List
        ) -> None:

        self.logger = logging.getLogger("all_in.Blastn")
        self.logger.debug("instantiating Blastn")
        self.settings = settings
        self.params = params
    
    def _execute(self, query: PathStr)-> List[BlastResultInfo]:
        """Execute blastn.

        Arguments:
            query(PathStr): path to the contig file.
        
        Returns:
            list of BlastResultInfo: results.
        
        Note:
            If the fasta file is like below:
            >OneSequence
            atgcatgc
            >anotherSequence
            atgcgcat

            The result is:
            [
                BlastResultInfo(OneSequence),
                BlastResultInfo(AnotherSequence)
            ] 
        """

        # Read query fasta
        with open(query, 'rt') as f:
            fasta = f.read().split('\n')[:-1]
            query_dict = {fasta[i]: fasta[i+1] for i in range(0,len(fasta),2)}
        self.logger.debug(query)
        self.logger.debug(fasta)
        self.logger.debug(query_dict)
        # Prepare result list
        result = list()

        # Search each sequence
        for name, seq in query_dict.items():
            # Work inside the context of tempfile.
            #with tempfile.TemporaryFile(mode="w+t") as tmp:
            with open("search_test.csv", "w+t") as tmp_out, tempfile.NamedTemporaryFile('w+t') as tmp_seq:
                # output sequence to tmp file
                tmp_seq.write(seq)
                tmp_seq.seek(0)
                self.logger.debug(f"query file: \n{tmp_seq.read()}")
                tmp_seq.seek(0)
                # Construct command
                command_line = [
                    "blastn",
                    "-query",
                    tmp_seq.name,
                    "-outfmt",
                    "10",
                    *self.params
                ]
                self.logger.debug("command_line: {}".format(command_line))

                self.logger.debug("execute {}".format(command_line))
                return_code = None
                try:
                    # Execute
                    # Input query seq(PIPE),
                    # Write result to temporary file(tmp),
                    # and get error from PIPE 
                    proc = subprocess.Popen(
                        command_line,
                        stdin=None,
                        stdout=tmp_out,
                        stderr=subprocess.PIPE
                    )
                    # Wait until search finished.
                    # We will get no stdout
                    _, err = proc.communicate()
                    return_code = proc.returncode
                    # Seek to file head
                    tmp_out.seek(0)

                except OSError as e:
                    self.logger.error(err)
                    self.logger.error("blastn returns {}".format(return_code))
                    self.logger.exception(e)
                    raise e
                else:
                    self.logger.debug(tmp_out.read())
                    tmp_out.seek(0)
                    # Convert csv to DataFrame
                    result_df = pd.read_csv(tmp_out, names=self.settings["blast_header"])
                    result.append(
                        BlastResultInfo(result_df,query,name,seq)
                    )
                    result_df.to_csv("1A01.csv")
        else:
            self.logger.debug("search end. result: {}".format(result))
        return result

    def _choose_highest_score(self, blastn_result: BlastResultInfo) -> BlastResultInfo:
        """Choose best result from blastn output.

        Arguments:
            blastn_result(BlastResultInfo): result of _execute.
        
        Returns:
            BlastResultInfo: best match.
        """
        # sort result by evalue
        blastn_result.result.sort_values(
            by=["evalue"], inplace=True)
        
        # Get highest evalue and extract result that have highest evalue
        highest_evalue = blastn_result.result["evalue"].min()
        highest_bitscore = blastn_result.result["bitscore"].max()
        self.logger.debug(f"eval={highest_evalue}, bitscore={highest_bitscore}")
        best_result: BlastResultInfo = BlastResultInfo(
            blastn_result.result.query(
                f"evalue == {highest_evalue} and bitscore == {highest_bitscore}"
            ),
            blastn_result.query_file,
            blastn_result.query_name,
            blastn_result.query_seq
        )

        return best_result

    def _get_result_intersection(self, results: List[BlastResultDF]) ->  Optional[BlastResultDF]:
        """Return intersection of results.

        Arguments:
            results(list of BlastResultDF): some results of blastn.
        
        Returns:
            BlastResultDF: intersection of results.
        
        Note:
            Returns can be one or None.
        """
        pass

    def _blast_search_single(self, query: PathStr) -> Optional[BlastResultInfo]:
        """Execute blast search and return best result.

        Arguments:
            query(PathStr): path to the query fasta file.
        
        Returns:
            BlastResultInfo: best match.
        
        Note:
            This is a private function.
            To execute search, use blast_search.
        """
        # Execute blastn
        self.logger.debug("_blast_search_single called")
        search_result: List[BlastResultInfo] = self._execute(query)

        # Extract highest score 
        if len(search_result) > 0:
            top_list: List[BlastResultInfo] = []
            for r in search_result:
                top_list.append(self._choose_highest_score(r))
            else:
                # Select df that has highest e-value among top_list 
                highest = min(
                    top_list,
                    key=(lambda x: x.result["evalue"].min())
                )
                res = highest
                self.logger.debug("top_list: {}".format(top_list))
                self.logger.debug("highest: {}".format(highest.result))
        else:
            # If there is no hit, return None
            res = None

        return res

    def _blast_search_multi(self, query_list: List[PathStr]) -> Optional[BlastResultInfo]:
        """Execute blast search and return best result.

        Arguments:
            query_list(list of PathStr): path to the query fasta file.
        
        Returns:
            BlastResultDF: best match.
        
        Note:
            This is a private function.
            To execute search, use blast_search.
        """
        pass

    def blast_search(self, query: Union[PathStr, List[PathStr]]) -> Optional[BlastResultInfo]:
        """Execute blast search and return best result as BlastResult(pd.Series)

        Arguments:
            query(PathStr or list of PathStr): path(s) to query fastq file(s).
        
        Returns:
            BlastResultDF: best match.
        """
        if isinstance(query, list):
            return self._blast_search_multi(query)
        elif isinstance(query, PathStr):
            return self._blast_search_single(query)
        else:
            self.logger.error("blast_search got {}: why?".format(type(query)))
            return None

if __name__ == "__main__":
    # test
    logging.basicConfig(level=logging.DEBUG, filename="blast.log", filemode="w")
    blast = Blastn(
        {"blast_header": [
            "qaccver",
            "saccver",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore"
        ]},
        [
        "-db",
        "USDA/USDA"
        ]
    )
    result = blast.blast_search(PathStr(sys.argv[1]))
    print(result)
    result.result.to_csv("result.csv")