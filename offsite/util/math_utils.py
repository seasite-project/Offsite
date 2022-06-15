"""@package util.math_utils
Util mathematical functions.

@author: Johannes Seiferth
"""

from typing import List, Optional, Tuple, Union

from numpy import asarray
from sympy import sympify, Symbol, solve


def eval_math_expr(expression: Union[str, float, int],
                   known_variables: Optional[List[Tuple[str, Union[str, float, int]]]] = None, cast_to=None):
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


def solve_equation(
        lhs: str, rhs: str, unknown: str, known_variables: Optional[List[Tuple[str, Union[str, float, int]]]] = None):
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


def remove_outliers(data: List[float]) -> List[float]:
    """Remove outliers from a list of points and return the filtered list.

    Parameters:
    -----------
    data: list of float
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
    observed_value: float
        Observed value.
    expected_value: float
        Expected value.

    Returns:
    --------
    float
        Percent deviation of the observed value from the expected value.
    """
    return (observed_value - expected_value) / expected_value * 100
