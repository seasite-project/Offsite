"""@package train.train_utils
Util functions used for training the tuning database.

@author: Johannes Seiferth
"""

from copy import deepcopy
from itertools import product
from sys import maxsize as sys_maxsize
from typing import Dict, List, Optional, Tuple

from pandas import DataFrame
from sqlalchemy.orm import Session
from sympy import simplify

import offsite.config
from offsite.config import Config
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel, KernelRecord
from offsite.train.impl_variant import ImplVariant
from offsite.util.sample_interval import SampleInterval, SampleType

IntervalRecordList = List[Tuple[SampleInterval, str]]


def fuse_equal_records(records: IntervalRecordList) -> IntervalRecordList:
    """
    Reduce total number of sample intervals by combining adjacent intervals. Adjacent intervals can be combined if they
    give the same result (e.g prediction, ...).

    Parameters:
    -----------
    records: list of tuple
        Results for a set of SampleInterval objects.

    Returns:
    --------
    list of tuple :
        Reduced set of the input SampleInterval objects.
    """
    if len(records) <= 1:
        return records
    records_cpy = deepcopy(records)
    # Fuse neighbouring intervals that give the same prediction.
    fused = list()
    interval = None
    result = -1
    for record in records_cpy:
        if record[0].type == SampleType.MODEL_INNER:
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


def fetch_kernel_runtime_prediction_data(db_session: Session, machine: int, compiler: int, method: int, ivp: int,
                                         cores: int, frequency: float, mode: Optional[str] = None,
                                         ode_size: Optional[int] = None) -> DataFrame:
    """
    Fetch all kernel runtime prediction data records available in the database that match the given configuration of
    machine configuration, frequency, cores, compiler, ODE method and IVP.

    If a ODE system size 'n' is passed only records for that fixed 'n' are returned.

    If a prediction mode is passed only records obtained using that mode are returned.

    Return these records sorted by the kernel id.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: int
        Used MachineState.
    compiler: int
        Used compiler.
    method: int
        Used ODE method.
    ivp: int
        Used IVP.
    cores: int
        Cores per NUMA domain of the used machine.
    frequency: float
        Clock frequency of the used machine.
    mode: str
        Application mode used to obtain data.
    ode_size: int
        Used fixed ODE system size.

    Returns:
    --------
    pandas.DataFrame
        Sorted kernel runtime prediction data records.
    """
    # Fetch data from database.
    data: DataFrame = KernelRecord.select(db_session, machine, compiler, method, ivp, cores, frequency, mode, ode_size)
    data.sort_values(['kernel', 'first', 'last'], inplace=True)
    return data


def deduce_available_impl_variants(db_session: Session, skeleton: ImplSkeleton, exclude_unchanged: bool,
                                   predicted_kernels: Optional[List[Kernel]] = list(),
                                   filter_yasksite_opt: bool = False) -> Tuple[List[ImplVariant], Dict[int, str]]:
    """
    Deduce all available implementation variants of a given ImplSkeleton object. Therefore, the ids of the skeleton's
    KernelTemplates are selected and these are in turn used to fetch the ids of each template's actual kernel variants.
    The cartesian product of the KernelTemplates' permutations (of possible kernel variants) determines the set of
    possible implementation variants.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    skeleton: ImplSkeleton
        Used ImplSkeleton object.
    predicted_kernels: List of Kernel
        List of kernels for which new predictions were obtained in this run.

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
    tid = -1
    for connect in skeleton.connected_templates:
        tid = connect[0].db_id
        # Select IDs of all kernel associated with kernel template 'tid'.
        kernel_ids[tid] = Kernel.select_kernels(db_session, tid)
    # We do not have to recompute the implementation predictions of a variant if the predictions of none of its
    # kernels changed! Thus, we ignore such variants if exclude_unchanged is set.
    if exclude_unchanged:
        no_changes = True
        only_one_changed = False
        changed = None
        for kids in kernel_ids.values():
            # Check for newly predicted kernels. Filter out those kernels for which no new predictions were obtained in
            # this run.
            tmp = [(kid[0], kid[1]) for kid in kids if kid[0] in predicted_kernels]
            # Both lists only have the same length if we obtained new predictions for all kernels. This means no
            # kernels could be filtered out.
            if len(tmp) != len(kids):
                # Otherwise we have to check if we filtered out kernels of more than one of the associated kernel
                # templates. We can skip variants only if just one kernel template changed.
                if no_changes:
                    changed = (tid, tmp)
                    only_one_changed = True
                    no_changes = False
                else:
                    only_one_changed = False
        # Return a empty list of variants since no kernels were predicted anew in this run. Thus, we do not need to
        # update the implementation variant predictions either!
        if len(predicted_kernels) == 0:
            return list(), dict()
        # Predictions of only one kernel template changed. Thus, we might be able to skip all implementation variants
        # that include kernels that were not predicted anew.
        if only_one_changed:
            assert changed is not None
            # At least one kernel of each kernel template is required to derive valid implementation variants.
            if len(changed[1]) > 0:
                kernel_ids[changed[0]] = changed[1]
    # Deduce all available implementation variant permutations from the given kernel ids.
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
    # Create all implementation variant objects.
    impl_variants: List[ImplVariant] = list()
    loaded_kernels: Dict[int, Kernel] = dict()
    for variant in variants:
        clust_lvl_comm, node_lvl_comm, loaded_kernels = ImplVariant.count_kernel_communication(
            db_session, variant, loaded_kernels)
        impl_variants.append(ImplVariant(skeleton.db_id, variant, clust_lvl_comm, node_lvl_comm))
    # Determine the number of times the particular templates are executed.
    executions = {kid[0]: connect[1] for connect in skeleton.connected_templates for kid in
                  kernel_ids[connect[0].db_id]}
    return impl_variants, executions


def deduce_impl_variant_sample_intervals(data: DataFrame) -> Dict[SampleInterval, Tuple[float, float]]:
    """Deduce the set of sample intervals of an implementation variant from its given kernel runtime prediction data.

    Parameters:
    -----------
    data: DataFrame
        Kernel runtime prediction data.

    Returns:
    --------
    dict of SampleInterval
        Dict of sample intervals.
    """
    config: Config = offsite.config.offsiteConfig
    # Copy prediction data to handle data in function.
    pdata = {kid: list() for kid in data['kernel'].unique()}
    for record in data.itertuples():
        pdata[record.kernel].append((record.first, record.last, record.prediction, record.weight))
    # Deduce samples intervals.
    impl_var_data = dict()
    # Fixed ODE system size.
    if config.args.ode_size:
        ode_size = config.args.ode_size
        impl_var_data[SampleInterval(ode_size, ode_size, None)] = dict()
        for kernel in pdata:
            impl_var_data[SampleInterval(ode_size, ode_size, None)][kernel] = (pdata[kernel][0][2], pdata[kernel][0][3])
        return impl_var_data
    end = -1
    while any(len(le) > 1 for le in pdata.values()):
        start = 1
        end = sys_maxsize
        # Determine lowest end value of current samples.
        low_sample = None
        for data in pdata.values():
            start = data[0][0]
            low_sample = data[0] if data[0][1] < end else low_sample
            end = min(data[0][1], end)
        # Add sample interval with this end value to 'impl_var_data'.
        impl_var_data[SampleInterval(start, end, None)] = dict()
        for kernel in pdata:
            impl_var_data[SampleInterval(start, end, None)][kernel] = (pdata[kernel][0][2], pdata[kernel][0][3])
        # Update first of all current samples.
        for data in pdata.values():
            data[0] = (end + 1, data[0][1], data[0][2], data[0][3])
        # Remove all current samples already represented in the list of significant samples.
        for data in pdata.values():
            if len(data) > 1 and data[0][0] > data[0][1]:
                data.pop(0)
    # Add interval with 'inf' end.
    end = max(end, 1)
    #
    impl_var_data[SampleInterval(end, sys_maxsize, None)] = dict()
    for kernel in pdata:
        impl_var_data[SampleInterval(end, sys_maxsize, None)][kernel] = (pdata[kernel][0][2], pdata[kernel][0][3])
    return impl_var_data
