"""filtering the result of blast.py
"""

import logging
import os
import sys
from typing import Dict, List, Optional

import pandas as pd
from tqdm.std import tqdm

if __name__ == "__main__":
    from my_types import *
else:
    from all_in_tools.my_types import *

def save_result(
    blast_result: Dict[str, Optional[BlastResultInfo]],
    out_csv_path: PathStr,
    filter_query: Optional[str] = None
    ) -> None:
    """Convert dict of BlastResultInfo into readable csv.

    Arguments:
        blast_result(Dict[str, Optional[BlastResultInfo]]): returns from blast.blast_*
        out_csv_path(PathStr): path to result.csv
        filter_query(Optional[str]): query for filtering result
    """
    logger = logging.getLogger("all_in.save_result")
    logger.debug("save_result called")
    result_df = pd.DataFrame(
        None, columns=[
            "plate", "cell", "candidate", "percent.ident.",
            "length", "evalue", "bitscore", "query_seq",
            "query_file", "from_intersection", "raw_count",
            "query_count", "other_candidates"
        ]
    )

    for cell, info in blast_result.items():
        # Continue if there is no result
        if info is None:
            continue

        # Get cell path and count raw sequences
        if info.intersection:
            cell_path = os.path.dirname(info.query_file.split(':')[0])
        else:
            cell_path = os.path.dirname(info.query_file)

        with open(os.path.join(cell_path, "tmp_R1.fastq"), "rt") as f:
            raw_count = len(f.read().split('\n'))//4
        
        # Count query sequences
        if info.intersection:
            query_count = 0
            for q in info.query_file.split(':'):
                with open(q, "rt") as f:
                    query_count += len(f.read().split('\n'))//2
        else:
            with open(info.query_file, "rt") as f:
                query_count = len(f.read().split('\n'))//2
        

        # Separate top and other candidate
        top_candidate = info.result.iloc[0,]
        other_candidate = info.result.iloc[1:,]

        # Construct row of result_df
        row_dict = {
            "plate": cell[0],
            "cell": cell[1:],
            "candidate": top_candidate["saccver"],
            "percent.ident.": top_candidate["pident"],
            "length": top_candidate["length"],
            "evalue": top_candidate["evalue"],
            "bitscore": top_candidate["bitscore"],
            "query_seq": top_candidate["query_seq"],
            "query_file": top_candidate["query_file"],
            "raw_count": raw_count,
            "query_count": query_count,
            "from_intersection": info.intersection,
            "other_candidates": list(other_candidate.loc[:,"saccver"])
        }
        row_series = pd.Series(row_dict)
        logger.debug(row_series)
        result_df = result_df.append(row_series, ignore_index=True)
    else:
        # Sort result
        result_df.sort_values(["plate", "cell"], inplace=True)

        # Filter result if filter is given
        result_df.query(filter_query, inplace=True) if filter_query is not None else None

        # write result
        if out_csv_path is None:
            out_csv_path = sys.stdout
        try:
            result_df.to_csv(out_csv_path, index=False)
        except OSError as e:
            logger.exception(e)
        else:
            logger.debug(f"write csv to {out_csv_path}")