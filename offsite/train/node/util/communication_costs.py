"""@package train.node.util.node_communication
Functions to train the tuning database with communication costs.

@author: Johannes Seiferth
"""

from typing import Dict, Optional

from pandas import DataFrame

from offsite.descriptions.ode import IVP, ODEMethod, ivp_grid_size, corrector_steps, stages
from offsite.util.math_utils import eval_math_expr


def compute_node_lvl_communication_costs(bench_data: DataFrame, comm_ops_node_lvl: Dict[str, str], max_cores: int,
                                         method: Optional[ODEMethod], ivp: Optional[IVP]) -> Dict[int, float]:
    """
    Compute the node-level communication costs when executing a given set of communication operations for a given
    configuration of ODE method and IVP.

    Parameters:
    -----------
    bench_data: pandas.DataFrame
        Available benchmark data sorted by number of cores and benchmark name.
    comm_ops_node_lvl: dict (key=str, value=str)
        Executed communication operations.
    method: ODEMethod
        Used ODE method.
    ivp: IVP
        Used IVP.

    Returns:
    dict(key: number of cores)
        Communication costs sorted by number of cores.
    """
    constants = list()
    if method is not None:
        constants.extend([corrector_steps(method), stages(method), ])
    if ivp is not None:
        constants.append(ivp_grid_size(ivp.gridSize))
    # Init empty dictionary in case we have no benchmark data.
    costs: Dict[int, float] = {cores: 0.0 for cores in range(1, max_cores)}
    # Compute for all number of cores...
    for _, row in bench_data.iterrows():
        operation = row['name']
        # Select measured benchmark result for this communication operation and multiply it with the number of
        # repetitions per iteration step of this operation.
        c = 0.0
        if operation in comm_ops_node_lvl:
            c += row['data'] * eval_math_expr(comm_ops_node_lvl[operation], constants, cast_to=float)
        costs[row['cores']] = c
    return costs
