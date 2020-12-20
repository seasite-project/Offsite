"""@package performance_model
Performance modeling functions.
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple

import attr

from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.evaluation.math_utils import eval_math_expr, ivp_grid_size, ivp_system_size


def compute_pmodel_kernel_pred(num_iterations: int, machine: Machine, ecm: float) -> float:
    """Compute the pmodel kernel prediction in cycles of a PModelKernel.

    The kernel prediction states the total number of cycles executed when executing a given number of iterations of
    that pmodel kernel.

    Parameters:
    -----------
    num_iterations: int
        Number of iterations executed by the PModelKernel considered.
    machine: Machine
        Machine the kernel is run on.
    ecm: float
        ECM result (in cycles per cache line) obtained by the kerncraft tool for the PModelKernel considered when
        executing 'num_iterations' iterations.

    Returns:
    --------
    float
        Kernel prediction in cycles of a PModelKernel.
    """
    constants = [('ecm', ecm), ('iterations', num_iterations), ('elements_cacheline', machine.elements_per_cacheline)]
    return eval_math_expr('ecm * iterations / elements_cacheline', constants)


def compute_kernel_runtime_pred(pmodel_predictions: List[str], frequency: float) -> float:
    """Compute the kernel runtime prediction in seconds of a Kernel.

    The kernel runtime prediction states the runtime in seconds when executing that kernel with a particular CPU
    frequency.

    Parameters:
    -----------
    pmodel_predictions: list of str
        PModel kernel predictions obtained for the pmodel kernels associated with this Kernel when running on a
        particular number of CPU cores.
    frequency: float
        CPU frequency the kernel is executed with.

    Returns:
    --------
    float
        Kernel runtime prediction in seconds.
    """
    # Construct equation.
    known_variables = [('f', frequency)]
    equation_str = ''
    for i, prediction in enumerate(pmodel_predictions):
        equation_str += '+ p{0}'.format(i)
        known_variables.extend([('p{}'.format(i), prediction)])
    equation_str = '({})/f'.format(equation_str)
    return eval_math_expr(equation_str, known_variables)


def compute_impl_variant_runtime_pred(associated_kernels: Tuple[int], kernel_runtime_preds: Dict[int, float],
                                      executions: Dict[int, float], communication_costs: float) -> float:
    """Compute the node-level runtime prediction in seconds of an implementation variant.

    The node-level runtime prediction states the runtime in seconds when executing a single time step of that
    implementation variant. The prediction is computed using the kernel runtime predictions of its associated kernels
    as well as the communication costs of this implementation variant.

    Parameters:
    -----------
    associated_kernels: tuple of int
        Ids of the Kernel objects associated with this implementation variant.
    kernel_runtime_preds: dict (key: Kernel object id)
        Kernel runtime predictions obtained for the kernels associated with this implementation variant.
    executions: dict (key: Kernel object id)
        Number of times the kernels associated with this implementation variant are run when executing a single time
        step.
    communication_costs: float
        Communication costs per step of this implementation variant.

    Returns:
    --------
    float
        Node-level runtime prediction in seconds.
    """
    # Compute the node-level runtime prediction.
    node_level_prediction = eval_math_expr(0.0)
    # Add up the single kernel runtime predictions.
    for kernel in associated_kernels:
        try:
            expr = '{} * {} * {}'.format(
                kernel_runtime_preds[kernel][0], kernel_runtime_preds[kernel][1], executions[kernel])
            node_level_prediction += eval_math_expr(expr)
        except KeyError:
            raise RuntimeError('Failed to compute impl variant runtime prediction!')
    # Add communication costs.
    node_level_prediction += communication_costs
    return node_level_prediction


class SamplePosition(Enum):
    """
    Defines if the sample lies in the border region of two adjacent intervals or not.

    - INNER
        Sample interval in the inner region of an interval.
    - BORDER
        Sample interval in the border region of two adjacent intervals.
    """
    INNER = 'INNER'
    BORDER = 'BORDER'


@attr.s(hash=True)
class SampleInterval:
    """Representation of a SampleInterval object.

    Attributes:
    -----------
    first: int
        First value included in the sample interval.
    last: int
        First value included in the sample interval.
    sample: int
        Actual value used to sample the interval.
    position: SamplePosition
        Type of the interval. Can be either 'INNER' or 'BORDER'.
    """
    first = attr.ib(type=int)
    last = attr.ib(type=int)
    sample = attr.ib(type=int, default=None)
    region = attr.ib(type=SamplePosition, default=SamplePosition.INNER)

    def __attrs_post_init__(self):
        self.first = int(self.first)
        self.last = int(self.last)
        if self.sample is not None:
            self.sample = int(self.sample)

    def __hash__(self):
        return hash((self.first, self.last))

    def __eq__(self, other):
        return (self.first, self.last) == (other.first, other.last)

    def __ne__(self, other):
        return not self == other

    def median(self, ivp: Optional[IVP] = None) -> int:
        """Return median value of this interval.

        If an IVP is passed not the exact median is returned but the nearest value that satisfies the constraints of
        the IVP regarding possible system sizes(e.g. square system size for Heat2D, ...).

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        int
            Median value of this interval.
        """
        median = int(round(self.first + (self.last - self.first) / 2))
        # Consider restraints of the grid size of the IVP constraining possible but also sampleable, median values.
        if ivp:
            # Round up to previous ...
            grid_size = int(eval_math_expr(ivp.gridSize, [ivp_system_size(median)]))
            median_low = int(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)]))
            if median_low < self.first or median_low > self.last:
                # Round up to next ...
                grid_size = int(round(eval_math_expr(ivp.gridSize, [ivp_system_size(median)])))
                median_high = int(round(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)])))
                if median_high < self.first or median_high > self.last:
                    median = None
                else:
                    median = median_high
            else:
                median = median_low
        return median
