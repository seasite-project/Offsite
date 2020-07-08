"""@package train_utils
Utility functions used for training the tuning database.
"""

from itertools import product
from sys import maxsize as sys_maxsize
from typing import Dict, List, Tuple

from sqlalchemy.orm import Session
from sympy import simplify

import offsite.config
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl_variant import ImplVariant
from offsite.descriptions.kernel_template import Kernel
from offsite.evaluation.math_utils import eval_math_expr, interpolate_polynomial, ivp_system_size
from offsite.evaluation.performance_model import SampleInterval, SamplePosition
from offsite.evaluation.records import KernelRecord

IntervalRecordList = List[Tuple[SampleInterval, str]]


def fuse_equal_records(records: IntervalRecordList) -> IntervalRecordList:
    """
    Reduce total number of sample intervals by combining adjacent intervals. Adjacent intervals can be combined if they
    give the same result (e.g prediction, ...).

    Parameters:
    -----------
    records : list of tuple
        Results for a set of SampleInterval objects.

    Returns:
    --------
    list of tuple :
        Reduced set of the input SampleInterval objects.
    """
    if len(records) <= 1:
        return records
    # Fuse neighbouring intervals that give the same prediction.
    fused = list()
    interval = None
    result = -1
    for record in records:
        if record[0].region == SamplePosition.INNER:
            if interval:
                fused.append((interval, result))
            interval = record[0]
            result = record[1]
        else:
            # Test if both intervals have the same prediction.
            if simplify(result - record[1]) == 0:
                # Increase upper bound of current interval.
                interval.last = record[0].last
            else:
                # Save interval.
                if interval:
                    fused.append((interval, result))
                # Update current interval and result.
                interval = record[0]
                result = record[1]
    # Save final interval.
    fused.append((interval, result))
    return fused


def interpolate_border_records(records: IntervalRecordList) -> IntervalRecordList:
    """Reduce total number of sample intervals by replacing border intervals with a single interpolated interval.

    Parameters:
    -----------
    records : list of tuple
        Results for a set of SampleInterval objects.

    Returns:
    --------
    list of tuple :
        Reduced set of the input SampleInterval objects.
    """
    if len(records) <= 1:
        return records
    # Interpolate border intervals.
    interpolated = list()
    border = list()
    for record in records:
        if record[0].region == SamplePosition.INNER:
            # Interpolate previous border region.
            if border:
                # Collect set of points that characterize the border region. Prediction is divided by system size 'n'
                # for the interpolation - later added again.
                # ... end previous inner region.
                x_coords = [interpolated[-1][0].last]
                y_coords = [eval_math_expr(interpolated[-1][1], [ivp_system_size(1)], cast_to=float)]
                # ... data points in border region.
                for interval in border:
                    x_coords.extend([interval[0].first, interval[0].last])
                    y_coords.extend([eval_math_expr(interval[1], [ivp_system_size(1)], cast_to=float),
                                     eval_math_expr(interval[1], [ivp_system_size(1)], cast_to=float)])
                # ... start following inner region.
                x_coords.append(record[0].first)
                y_coords.append(eval_math_expr(record[1], [ivp_system_size(1)], cast_to=float))
                # Interpolate polynomial.
                polynomial = interpolate_polynomial(x_coords, y_coords, 2)
                # Multiply with system size 'n' again.
                polynomial = '({}) * n'.format(polynomial)
                interval = SampleInterval(x_coords[0] + 1, x_coords[-1] - 1, region=SamplePosition.BORDER)
                # Save interpolated border region interval.
                interpolated.append((interval, polynomial))
                # Reset border sample
                border = list()
            # Save current inner region interval.
            interpolated.append(record)
        else:
            # Add to current border region.
            border.append(record)
    # Last sample should be an inner region for the main memory level. Hence, the list of border samples should be empty
    # at this point.
    assert not border
    return interpolated


def reduce_records(records: IntervalRecordList) -> IntervalRecordList:
    """Reduce total number of sample intervals by combining adjacent intervals and interpolating border intervals.

    Parameters:
    -----------
    records : list of tuple
        Results for a set of SampleInterval objects.

    Returns:
    --------
    list of tuple
        Reduced set of the input SampleInterval objects.
    """
    config = offsite.config.offsiteConfig
    # Fuse neighbouring intervals that give the same prediction.
    records = fuse_equal_records(records)
    # Interpolate border intervals.
    if config.args.ws_interpolate:
        records = interpolate_border_records(records)
    return records


def fetch_and_sort_kernel_runtime_prediction_data(
        db_session: Session, machine: int, compiler: int, method: int, ivp: int, cores: int, frequency: float,
        mode: str, ode_size: int = None) -> Dict[int, KernelRecord]:
    """
    Fetch all kernel runtime prediction data records available in the database that match the given configuration of
    machine, frequency, cores, compiler, ODE method and IVP.

    If a ODE system size 'n' is passed only records for that fixed 'n' are returned.

    Return these records sorted by the kernel id.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : int
        Used machine.
    compiler : int
        Used compiler.
    method : int
        Used ODE method.
    ivp : int
        Used IVP.
    cores : int
        Cores per NUMA domain of the used machine.
    frequency : float
        Clock frequency of the used machine.
    mode : str
        Application mode used to obtain data.
    ode_size : int
        Used fixed ODE system size.

    Returns:
    --------
    dict (key: kernel id, value: KernelRecord)
        Sorted kernel runtime prediction data records.
    """
    # Fetch data from database.
    data = KernelRecord.select(db_session, machine, compiler, method, ivp, cores, frequency, mode, ode_size)
    # Sort data.
    sorted_data = dict()
    for record in data:
        kernel = record.kernel
        # Sort data by kernel:
        if kernel in sorted_data:
            sorted_data[kernel].append(record)
        else:
            sorted_data[kernel] = [record]
    return sorted_data


def deduce_available_impl_variants(db_session: Session, skeleton: ImplSkeleton,
                                   filter_yasksite_opt: bool = False) -> Tuple[List[ImplVariant], Dict[int, str]]:
    """
    Deduce all available implementation variants of a given ImplSkeleton object. Therefore, the ids of the skeleton's
    KernelTemplates are selected and these are in turn used to fetch the ids of each template's actual kernel variants.
    The cartesian product of the KernelTemplates' permutations (of possible kernel variants) determines the set of
    possible implementation variants.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    skeleton : ImplSkeleton
        Used ImplSkeleton object.

    Returns:
    --------
    list of ImplVariant
        Set of possible implementation variants.
    dict (key: kernel id, value: str)
        Number of times the kernel variants, which can be part of the ImplSkeleton's variants, are executed when
        computing a single iteration step.
    """
    # Fetch the kernel variant ids from the database.
    kernel_ids = dict()
    for connect in skeleton.connected_templates:
        tid = connect[0].db_id
        kernel_ids[tid] = Kernel.select_kernels(db_session, tid)
    # Deduce all implementation variant permutations.
    variants = product(*kernel_ids.values())
    # Filter out all implementation variants whose kernels do not use the same optimization parameters.
    if filter_yasksite_opt:
        filtered_variants = list()
        for variant in variants:
            base_blocking = variant[0][1]['ys_blocking']
            base_folding = variant[0][1]['ys_folding']
            select_variant = True
            filtered_kernel_ids = list()
            for kernel in variant:
                # Filter based on optimization parameters.
                if (kernel[1]['ys_blocking'] != base_blocking) or (kernel[1]['ys_folding'] != base_folding):
                    select_variant = False
                    break
                filtered_kernel_ids.append(kernel[0])
            if select_variant:
                filtered_variants.append(filtered_kernel_ids)
        variants = filtered_variants
    else:
        # Convert to lists of kernel IDs associated with a particular variant.
        variants = [[kernel[0] for kernel in variant] for variant in variants]
    # Deduce all implementation variant objects.
    variants = [ImplVariant(skeleton.db_id, variant) for variant in variants]
    # Determine the number of times the particular templates are executed.
    executions = {kid[0]: connect[1] for connect in skeleton.connected_templates for kid in
                  kernel_ids[connect[0].db_id]}
    return variants, executions


def deduce_impl_variant_sample_intervals(
        data: Dict[int, List[KernelRecord]]) -> Dict[SampleInterval, Tuple[float, float]]:
    """Deduce the set of sample intervals of an implementation variant from its given kernel runtime prediction data.

    Parameters:
    -----------
    data : dict of list of KernelRecord
        Kernel runtime prediction data.

    Returns:
    --------
    dict of SampleInterval
        Dict of sample intervals.
    """
    config = offsite.config.offsiteConfig
    # Copy prediction data to handle data in function.
    # Deepcopy not possible here!
    # Still there should be a better solution!
    pred_data = dict()
    for kid, kdata in data.items():
        pred_data[kid] = list()
        for record in kdata:
            pred_data[kid].append((record.first, record.last, record.prediction, record.weight))
    # Deduce samples intervals.
    impl_var_data = dict()
    # Fixed ODE system size.
    if config.args.ode_size:
        ode_size = config.args.ode_size
        impl_var_data[SampleInterval(ode_size, ode_size)] = dict()
        for kernel in pred_data:
            impl_var_data[SampleInterval(ode_size, ode_size)][kernel] = (pred_data[kernel][0][2],
                                                                         pred_data[kernel][0][3])
        return impl_var_data
    end = -1
    while any(len(l) > 1 for l in pred_data.values()):
        start = 1
        end = sys_maxsize
        # Determine lowest end value of current samples.
        low_sample = None
        for pdata in pred_data.values():
            start = pdata[0][0]
            low_sample = pdata[0] if pdata[0][1] < end else low_sample
            end = min(pdata[0][1], end)
        # Add sample interval with this end value to 'impl_var_data'.
        impl_var_data[SampleInterval(start, end)] = dict()
        for kernel in pred_data:
            impl_var_data[SampleInterval(start, end)][kernel] = (pred_data[kernel][0][2], pred_data[kernel][0][3])
        # Update first of all current samples.
        for pdata in pred_data.values():
            pdata[0] = (end + 1, pdata[0][1], pdata[0][2], pdata[0][3])
        # Remove all current samples already represented in the list of significant samples.
        for pdata in pred_data.values():
            if len(pdata) > 1 and pdata[0][0] > pdata[0][1]:
                pdata.pop(0)
    # Add interval with 'inf' end.
    impl_var_data[SampleInterval(end, sys_maxsize)] = dict()
    for kernel in pred_data:
        impl_var_data[SampleInterval(end, sys_maxsize)][kernel] = (pred_data[kernel][0][2], pred_data[kernel][0][3])
    return impl_var_data
