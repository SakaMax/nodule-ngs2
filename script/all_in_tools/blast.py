"""Blastn wrapper for all_in.py
"""
from functools import reduce
from glob import glob
import logging
import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Set, Union
import sys

import pandas as pd
from tqdm.std import tqdm

if __name__ == "__main__":
    from my_types import *
else:
    from all_in_tools.my_types import *

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
                    # Add origin of result to df
                    result_df["query_file"] = query
                    result_df["query_name"] = name
                    result_df["query_seq"] = seq
                    # Append result of this execution into return list
                    result.append(
                        BlastResultInfo(result_df,query,name,seq)
                    )
                    # check result (for debug)
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

    def _get_result_intersection(self, results: List[BlastResultInfo]) ->  Optional[BlastResultInfo]:
        """Return intersection of results.

        Arguments:
            results(list of BlastResultDF): some results of blastn.
        
        Returns:
            BlastResultDF: intersection of results.
        
        Note:
            Returns can be one or None.
        """
        self.logger.debug("intersection called")
        # Extract DataFrames and  strain names
        # Each set contain names from one of the result
        dfs: List[pd.DataFrame] = [r.result for r in results]
        names: List[Set[str]] = [set(d["saccver"]) for d in dfs]

        # Create union for further convolution process
        union_names = set()
        for s in names:
            union_names = union_names | s

        # Get intersection by convolution
        # f_and = lambda x,y: x & y
        # intersection: List[str] = list(reduce(f_and, names, union_names))
        def f_and(x,y):
            inter = x & y
            self.logger.debug(inter)
            return inter
        self.logger.debug(f"initial (union={len(union_names)}) : {union_names}")
        intersection: List[str] = list(reduce(f_and, names, union_names))
        self.logger.debug(f"final (intersection={len(intersection)}): {intersection}")

        if len(intersection) > 0:
            # Prepare final result DataFrame
            final_cols = self.settings["blast_header"] + ["query_file", "query_name", "query_seq"]
            final_df = pd.DataFrame(index=[], columns=final_cols)

            # Append appropriate record("saccver" of record <= intersection ) into final_df
            # This is maybe a kind of convolution... 
            for d in dfs:
                final_df = pd.concat(
                    [final_df, d[d["saccver"].isin(intersection)]]
                )
            
            # Remove duplicated records
            # Save one has highest evalue
            final_df.sort_values(by=["evalue"], inplace=True)
            final_df.to_csv("final_before_drop.csv")
            final_df.drop_duplicates(subset="saccver", keep="first", inplace=True)

            # set query_file, query_name, and query_seq by convined name
            final_query_file = ":".join(list(final_df["query_file"].unique()))
            final_query_name = ":".join(list(final_df["query_name"].unique()))
            final_query_seq  = ":".join(list(final_df["query_seq"].unique()))

            # return BlastResultInfo
            return BlastResultInfo(
                final_df,
                final_query_file,
                final_query_name,
                final_query_seq,
                True
            )
        else:
            # No result
            return None

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
                # skip empty dataframe
                if r.result.empty:
                    continue
                top_list.append(self._choose_highest_score(r))
            else:
                if top_list: # if there is more than zero result(s)
                    # Select df that has highest e-value among top_list 
                    highest = min(
                        top_list,
                        key=(lambda x: x.result["evalue"].min())
                    )
                    res = highest
                    self.logger.debug("top_list: {}".format(top_list))
                    self.logger.debug("highest: {}".format(highest.result))
                else: # If there is no valid result
                    res = None
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
        # Call _blast_search_single for each file,
        # then get intersection

        self.logger.debug("_blast_search_multi called.")
        self.logger.debug(query_list)

        # List for result
        individual_result: List[BlastResultInfo] = list()

        for query in query_list:
            individual_result.append(
                self._blast_search_single(query)
            )
        else:
            # Remove None
            filter_none = lambda x: True if x is not None else False
            individual_result = list(filter(filter_none, individual_result))

        # Extract result from individual_result        
        final_result = None
        if len(individual_result) == 0:
            # There is no answer
            pass
        elif len(individual_result) == 1:
            # Exactly one result: return this
            final_result = individual_result[0]
        else:
            # There is more than one result
            # Get intersection
            final_result = self._get_result_intersection(individual_result)
        
        return final_result

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

def blast_all(
    destination: PathStr,
    settings: Dict
    ) -> Dict[str, Union[BlastResultInfo, None]]:
    """Blastn search, using all_contigs.fasta
    
    Arguments:
        destination(PathStr): path to the data folder
        settings(dict): settings from yaml

    Returns:
        dict[str: Optional[BlastResultInfo]]: cell_name: result of blast search
    """
    logger = logging.getLogger("all_in.blast_all")
    logger.debug("assemble_all called.")

    # Initialize blast

    blast = Blastn(
        settings=settings,
        params=settings["blastn"]
    )

    # Initialize result dict
    result_dict = dict()

    # Scan cells directory
    cells_dir = os.path.join(destination, "cells")
    with os.scandir(cells_dir) as it:
        for entry in tqdm(list(it), desc="blast search (all)"):
            if entry.is_dir():
                # Get cell name and the fasta file of final contigs
                cell_name = entry.path.split('/')[-1]
                contigs_path = os.path.join(cells_dir, cell_name, "all_contigs.fasta")

                # Execute search
                res = blast.blast_search(contigs_path)

                # Update result
                result_dict[cell_name] = res

    return result_dict

def blast_individual(
    destination: PathStr,
    settings: Dict
    ) -> Dict[str, Union[BlastResultInfo, None]]:
    """Blastn search, using *_ind_contigs.fasta
    
    Arguments:
        destination(PathStr): path to the data folder
        settings(dict): settings from yaml

    Returns:
        dict[str: Optional[BlastResultInfo]]: cell_name: result of blast search
    """
    logger = logging.getLogger("all_in.blast_individual")
    logger.debug("assemble_individual called.")

    # Initialize blast

    blast = Blastn(
        settings=settings,
        params=settings["blastn"]
    )

    # Initialize result dict
    result_dict = dict()

    # Scan cells directory
    cells_dir = os.path.join(destination, "cells")
    with os.scandir(cells_dir) as it:
        for entry in tqdm(list(it), desc="blast search (individual & intersection)"):
            if entry.is_dir():
                # Get cell name and the fasta file of final contigs
                cell_name = entry.path.split('/')[-1]
                contigs_path_list = glob(
                    os.path.join(cells_dir, cell_name, "*_ind_contigs.fasta")
                )

                # Execute blastn
                res = blast.blast_search(contigs_path_list)

                # Update result
                result_dict[cell_name] = res
    
    return result_dict

if __name__ == "__main__":
    # test
    logging.basicConfig(level=logging.DEBUG, filename="blast.log", filemode="w")
    mock_settings = {
        "blast_header": [
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
        ],
        "blastn": [
            "-db",
            "USDA/USDA"
            ]
        }

    b = Blastn(mock_settings, mock_settings["blastn"])
    print(b.blast_search("4C09.fasta"))

    #result_all = blast_all(sys.argv[1], mock_settings)
    #result_ind = blast_individual(sys.argv[1], mock_settings)
    #print(result_all)
    #print("--------")
    #print(result_ind)
    #result_ind["a_cell"].result.to_csv("result.csv")