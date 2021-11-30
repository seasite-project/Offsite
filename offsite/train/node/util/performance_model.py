"""@package train.node.util.performance_model
Performance modeling functions.
"""

from typing import Dict, List

from offsite.util.math_utils import eval_math_expr


def compute_pmodel_kernel_pred(num_iterations: int, ecm_pred: float, elements_per_cl: int) -> str:
    """Compute the pmodel kernel prediction (in cycles) of a pmodel kernel.

    The pmodel kernel prediction corresponds to the number of cycles required to execute a particular pmodel kernel
    for a given number of iterations and cacheline size.

    Parameters:
    -----------
    num_iterations: int
        Number of iterations executed by the pmodel kernel considered.
    ecm_pred: float
        ECM result (in cycles per cache line) obtained by the kerncraft tool for the PModelKernel considered when
        executing 'num_iterations' iterations.
    elements_per_cl: int
        xxx

    Returns:
    --------
    str
        Kernel prediction (in cycles) of a pmodel kernel.
    """
    known_vars = [('ecm', ecm_pred), ('iters', num_iterations), ('elements_cl', elements_per_cl)]
    pmodel_pred = eval_math_expr('ecm * iters / elements_cl', known_vars)
    return pmodel_pred


def compute_kernel_pred(pmodel_preds: List[str], frequency: float) -> str:
    """Compute the kernel runtime prediction (in seconds) of a kernel.

    The prediction depends on:
        (1) the pmodel kernel predictions of the associated pmodel kernels.
        (2) the CPU frequency used.

    The kernel runtime prediction corresponds to the runtime (in seconds) required to execute a particular
    kernel with a given CPU frequency.

    Parameters:
    -----------
    pmodel_preds: list of str
        PModel kernel predictions of the pmodel kernels associated with this kernel obtained on a particular
        machine and a particular number of CPU cores.
    frequency: float
        CPU frequency the kernel is executed with.

    Returns:
    --------
    str
        Kernel runtime prediction (in seconds) of a kernel.
    """
    known_vars = [('f', frequency)]
    # Construct equation...
    eqn_str = ''
    # ... sum up all pmodel predictions
    for idx, pred in enumerate(pmodel_preds):
        pi = 'p{}'.format(idx)
        eqn_str += '+ {}'.format(pi)
        known_vars.extend([(pi, pred)])
    # ... and divide through the frequency.
    eqn_str = '({})/f'.format(eqn_str)
    # Evaluate equation.
    kernel_pred = eval_math_expr(eqn_str, known_vars)
    return kernel_pred


FloatDict = Dict[int, float]


def compute_impl_variant_pred(kernel_preds: FloatDict, kernel_execs: FloatDict, comm_costs: float) -> str:
    """Compute the node-level runtime prediction (in seconds) of an implementation variant.

    The node-level runtime prediction corresponds to the runtime (in seconds) required to execute a particular
    implementation variant.

    The prediction depends on:
        (1) the kernel predictions of the associated kernels.
        (2) the communication costs of the implementation variant.

    Parameters:
    -----------
    kernel_preds: dict (key: Kernel object id)
        Kernel runtime predictions of the kernels associated with this implementation variant.
    kernel_execs: dict (key: Kernel object id)
        Number of times each associated kernel is executed when running this implementation variant.
    comm_costs: float
        Communication costs of this this implementation variant.

    Returns:
    --------
    str
        Node-level runtime prediction (in seconds) of an implementation variant.
    """
    # Sanity check.
    try:
        assert kernel_preds.keys() == kernel_execs.keys()
    except AssertionError:
        raise RuntimeError('Failed to compute the implementation variant runtime prediction!\nError: Keys should '
                           'match:\n* Kernel predictions: {}\n* Executions: {}'.format(kernel_preds, kernel_execs))
    # Compute the node-level runtime prediction as the sum of ...
    # ... the communication costs
    impl_pred: str = eval_math_expr(comm_costs)
    # ... and the single kernel runtime predictions.
    for kid, kpred in kernel_preds.items():
        expr: str = '{} * {} * {}'.format(kpred[0], kpred[1], kernel_execs[kid])
        impl_pred += eval_math_expr(expr)
    impl_pred = eval_math_expr(impl_pred)
    return impl_pred
