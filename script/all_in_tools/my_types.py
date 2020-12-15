from typing import NamedTuple, Callable

import pandas as pd

# Type alias
FASTA = str
PathStr = str
BlastResultDF = pd.DataFrame

# Tuple
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
    intersection: bool = False

class WorkflowFunctionInfo(NamedTuple):
    """The set of infomation AllIn.workflow_generator send.

    Attributes:
        function(callable): work in this step
        name(str): name of this step
        step(int): which step?
    """
    function: Callable
    name: str
    step: int