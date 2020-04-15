"""@package math_utils
Util mathematical functions.
"""

from typing import List, Tuple

from numpy import polyfit, asarray
from sympy import sympify, Symbol, solve

from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod


def corrector_steps(method: ODEMethod) -> Tuple[str, str]:
    """Return corrector steps variable.

    Parameters:
    -----------
    method : ODEMethod
        Used ODE method.

    Returns:
    --------
    tuple(str, str)
        Corrector step variable.
    """
    return 'm', method.correctorSteps


def stages(method: ODEMethod) -> Tuple[str, str]:
    """Return stages variable.

    Parameters:
    -----------
    method : ODEMethod
        Used ODE method.

    Returns:
    --------
    tuple(str, str)
        Stage variable.
    """
    return 's', method.stages


def ivp_system_size(system_size: str) -> Tuple[str, str]:
    """Return IVP system size variable.

    Parameters:
    -----------
    system_size : str
        System size of the solved IVP.

    Returns:
    --------
    tuple(str, str)
        IVP system size variable.
    """
    return 'n', system_size


def ivp_grid_size(grid_size) -> Tuple[str, str]:
    """Return IVP grid size variable.

    Parameters:
    -----------
    grid_size : str
        Grid size of the solved IVP.

    Returns:
    --------
    tuple(str, str)
        IVP grid size variable.
    """
    return 'g', grid_size


def cacheline_elements(machine: Machine) -> Tuple[str, str]:
    """Return number of elements per cache line variable.

    Parameters:
    -----------
    machine : Machine
        Used machine.

    Returns:
    --------
    tuple(str, str)
        Elements per cache line variable.
    """
    return 'c', machine.elements_per_cacheline


def eval_math_expr(expression: str, known_variables: List[Tuple[str, float]] = None, cast_to=None):
    """Evaluate mathematical expression string and return its value.

    Parameters:
    -----------
    expression: str
        Mathematical expression to be evaluated.
    known_variables: list of tuple(var, val)
        List of known variables of this expression.
    cast_to datatype
        Try to cast the evaluated mathematical expression to the passed data type. Might raise an error if the
        expression can not be casted to the passed type does not support conversions.

    Returns:
    --------
    sympy number/expression
        Value of the evaluated mathematical expression.
    """
    # Create sympy object.
    expression = sympify(expression)
    # Substitute known variables in equation.
    if known_variables:
        # Convert known_variables to sympy objects.
        known_variables_sympy = list()
        for var in known_variables:
            known_variables_sympy.append((var[0], sympify(var[1])))
        expression = expression.subs(known_variables_sympy)
    evaluated = expression.evalf()
    if cast_to:
        evaluated = cast_to(evaluated)
    return evaluated


def solve_equation(lhs: str, rhs: str, unknown: str, known_variables: List[Tuple[str, float]] = None):
    """Solve mathematical equation given as strings and return its solution.

    Parameters:
    -----------
    lhs: str
        Left-hand side of this equation.
    rhs: str
        Right-hand side of this equation.
    unknown: str
        Unknown of the solution.
    known_variables: list of tuple(var, val)
        List of known variables of this equation.

    Returns:
    --------
    list
        Solution(s) of this equation.
    """
    # Create sympy objects.
    equation = sympify('Eq({},{})'.format(lhs, rhs))
    unknown = Symbol(unknown)
    # Substitute known variables in equation.
    if known_variables:
        # Convert known variables to sympy objects.
        known_variables_sympy = list()
        for var in known_variables:
            known_variables_sympy.append((var[0], sympify(var[1])))
        equation = equation.subs(known_variables_sympy)
    # Solve equation.
    return solve(equation, unknown)


def interpolate_polynomial(x_coordinates: List[float], y_coordinates: List[float], degree: int) -> str:
    """Interpolate polynomial for a set of points and return the corresponding interpolation polynomial.

    Parameters:
    -----------
    x_coordinates : list of float
        Known x-coordinates. Must be sorted in increasing order.
    y_coordinates : list of float
        Known y-coordinates.
    degree : int
        Degree of the fitting polynomial.

    Returns:
    --------
    str
        Interpolation polynomial.
    """
    assert degree > 1
    # Interpolate the set of data points.
    coefficients = polyfit(x_coordinates, y_coordinates, degree)
    # Create polynomial expression.
    polynomial = ''
    for idx, coefficient in enumerate(coefficients):
        polynomial += ' + {} * x^{}'.format(coefficient, degree - idx)
    # Simplify polynomial expression if possible.
    polynomial = eval_math_expr(polynomial)
    return polynomial


def remove_outliers(data: List[float]) -> List[float]:
    """Remove outliers from a list of points and return the filtered list.

    Parameters:
    -----------
    data : list of float
        Data points with outliers.

    Returns:
    --------
    list of float
        Data minus outliers.
    """
    # Nothing do if equal data points.
    if all(x == data[0] for x in data):
        return data
    # Convert data to an numpy array.
    array = asarray(data)
    # Compute mean and standard deviation of array.
    mean = array.mean()
    std = array.std()
    # Remove outliers.
    data_filtered = [x for x in array if x > mean - 2 * std]
    data_filtered = [x for x in data_filtered if x < mean + 2 * std]
    return data_filtered


def percent_deviation(observed_value: float, expected_value: float) -> float:
    """Compute percent deviation of an observed value from the expected value.

    Parameters:
    -----------
    observed_value : float
        Observed value.
    expected_value : float
        Expected value.

    Returns:
    --------
    float
        Percent deviation of the observed value from the expected value.
    """
    return (observed_value - expected_value) / expected_value * 100
