"""@package util.working_sets.direct_formula_defintion
Definition of class DirectFormulaDefinition.

@author: Hana Shatri
"""

from re import fullmatch, sub, split

from numpy import array_equal
from pandas import DataFrame
from sympy import Symbol, simplify, sympify
from sympy.utilities.iterables import variations

from offsite.codegen.code_dsl.code_node import CodeNodeType


class DirectFormulaDefinition:

    def __init__(self, code_tree_in_df):
        self.code_tree_in_df = code_tree_in_df

    def unique_references_count(self, difficult_cases=None) -> str:
        """ Returns directly the formula.
        """
        if difficult_cases is None:
            difficult_cases = []
        computations = self.code_tree_in_df.loc[(self.code_tree_in_df["type"] == CodeNodeType.COMPUTATION)]
        references_list = []

        # formula
        count_total_references = 0
        for index, computation_node in computations.iterrows():
            for reference in computation_node["references"]:
                # surrounding_loops_df = reference["surrounding_loops"]
                if not self._check_if_reference_exists(reference, references_list) and \
                        not difficult_cases.__contains__(reference["reference_name"]):
                    # check dimension of the array, multiply accesses based on type of expression
                    count_current_refs = self._accesses_though_all_dimensions(reference)
                    references_list.append(reference)
                    count_total_references += count_current_refs

        return count_total_references

    # returns the formula
    def non_unique_simple_references_count(self) -> str:
        computations = self.code_tree_in_df.loc[(self.code_tree_in_df["type"] == CodeNodeType.COMPUTATION)]

        references_list = []
        # unify unique references
        for index, computation_node in computations.iterrows():
            for reference in computation_node["references"]:
                # if not references_list.__contains__(ref):
                if not self._check_if_reference_exists(reference, references_list):
                    references_list.append(reference)

        # group references by their name
        references_groups = DataFrame(references_list).groupby("reference_name")

        # formula
        count_total_references = 0
        # for each group of references with the same name
        for name, group in references_groups:

            if len(group) == 1:
                # unique references are easily calculated using the boundary variables.
                count_current_refs = 1
                # if the reference is multidimensional, then multiply the results.
                surrounding_loops_df = group.iloc[0]["surrounding_loops"]
                for arg_exp in group.iloc[0]["arg_expressions"]:
                    count_current_refs *= self._accesses_through_one_dimension(arg_exp, surrounding_loops_df)
                # add the results of one reference to the total.
                count_total_references += count_current_refs
            else:
                count_current_refs = self._equivalent_references_with_different_vars_same_bounds(group)

                if count_current_refs is None:
                    # non-unique references -> a group of more than one reference with the same name,
                    # but different arg expressions
                    argument_expressions = []
                    # for every reference, do the following
                    # if we have: ref[k-1][k+1], ref[k+2][k+3]
                    # then we have
                    # [{1: k-1, 2: k+1, "surrounding_loops": df},
                    #  {1: k+2, 2: k+3, "surrounding_loops": df}]
                    i = 0
                    for index, reference in group.iterrows():
                        i = 0
                        arg_dict = {}
                        for exp in reference["arg_expressions"]:
                            # ref[k+1][k-3][k] -> {1: k+1, 2: k-3, 3: k}
                            arg_dict[i] = exp
                            i += 1
                        arg_dict["surrounding_loops"] = reference["surrounding_loops"]
                        argument_expressions.append(arg_dict)

                    df_argument_expressions = DataFrame(argument_expressions)
                    count_current_refs = 1
                    # find biggest and smallest +- val for the same positional argument
                    # format "var +- value"
                    for columns in range(0, i):
                        # unique cases
                        if len(df_argument_expressions.groupby(columns)) == 1:
                            surrounding_loops_df = df_argument_expressions["surrounding_loops"].values[0]
                            count_current_refs *= self._accesses_through_one_dimension(
                                df_argument_expressions[columns].values[0], surrounding_loops_df)

                        else:
                            # the condition is for the variables of the same dimension to go up to the same value,
                            # so it does not matter which var value is used, because all are the same
                            # minimum is created based on the arg expression structure
                            # here we treat: k+-num, num*k+-num, k*num+-num, k/num+-num, k%num+-num
                            # first case: k+-num, k%num+-num
                            # second case: num*k+-num, k*num+-num
                            # third case: k/num+-num
                            var_boundary = None
                            list_of_expressions = []
                            exp_type = None
                            for index, elements in df_argument_expressions[[columns, "surrounding_loops"]].iterrows():
                                # i, 0-> n;  i + 10, i - 10,  = n + | 20 - (-20) | = min(n + 40, 2*n)
                                # n=5,   i+ 0, i+6 = min(n+6, 2*n) = min(11, 10) = 10
                                # 2*i + 4, 3*i +5 # todo handle gcd cases
                                el = str(simplify(elements[columns])).replace(" ", "")
                                surrounding_loops_df = elements["surrounding_loops"]
                                var, value_at_var, int_val2, sign, exp_type = _trivial_bound(el)
                                var_boundary = self._accesses_through_one_dimension(el, surrounding_loops_df)
                                list_of_expressions.append({"var": var,
                                                            "value_at_var": value_at_var,
                                                            "value": float(sign + str(int_val2)),
                                                            "type": exp_type})

                            # df_of_expressions = pd.DataFrame(list_of_expressions)
                            if exp_type == "simple" or exp_type == "mult_vars":
                                minimum_expression = self._min_expression(var_boundary, DataFrame(list_of_expressions))
                            else:
                                minimum_expression = self._different_cases_of_multiplication(var_boundary,
                                                                                             list_of_expressions)
                            try:
                                count_current_refs *= eval(minimum_expression)
                            except:
                                count_current_refs *= minimum_expression

                if count_current_refs != 1:
                    # add the results of one reference to the total.
                    count_total_references += count_current_refs
        count_total_references = simplify(count_total_references)
        return str(count_total_references)

    def _equivalent_references_with_different_vars_same_bounds(self, same_name_references_group: DataFrame):
        # replace all argument variables of the first reference with boundary values
        first_reference = same_name_references_group.iloc[0]
        replaced_expressions = []
        for exp in first_reference["arg_expressions"]:
            variable, _, _, _, _ = _trivial_bound(exp)
            replaced_exp = (str(exp)).replace(variable, str(
                first_reference["surrounding_loops"][(first_reference["surrounding_loops"]["var"] ==
                                                      variable)]["boundary"].iloc[0]))
            replaced_expressions.append(replaced_exp)

        for index, reference in same_name_references_group.iterrows():
            temp_replaced_expressions = []
            for exp in reference["arg_expressions"]:
                variable, _, _, _, _ = _trivial_bound(exp)
                replaced_exp = str(exp).replace(variable, str(
                    reference["surrounding_loops"][(reference["surrounding_loops"]["var"] ==
                                                    variable)]["boundary"].iloc[0]))
                temp_replaced_expressions.append(replaced_exp)

            for i in range(len(replaced_expressions)):
                if simplify(sympify(replaced_expressions[i]) - sympify(temp_replaced_expressions[i])) \
                        != 0:
                    return None

        return self._accesses_though_all_dimensions(first_reference)

    @staticmethod
    def _accesses_through_one_dimension(arg_exp, surrounding_loops_df):
        # variable, int multiplied with variable, int after +/- sign, sign
        if arg_exp.isdigit():
            return 1
        var, value_at_var, _, _, exp_type = _trivial_bound(arg_exp)

        if exp_type == "simple":
            if value_at_var is None:
                return surrounding_loops_df[surrounding_loops_df["var"] == var]["boundary"].values[0]
            elif value_at_var is not None and var is not None:
                return surrounding_loops_df[surrounding_loops_df["var"] == var]["boundary"].values[0] * value_at_var
            else:
                # modulo cases
                return value_at_var
        elif exp_type == "mult":
            if value_at_var < 0:
                return surrounding_loops_df[surrounding_loops_df["var"] == var]["boundary"].values[0] * value_at_var
            else:
                return surrounding_loops_df[surrounding_loops_df["var"] == var]["boundary"].values[0]
        else:
            result = 0
            for v in var:
                if not surrounding_loops_df[surrounding_loops_df["var"] == v].empty:
                    result += surrounding_loops_df[surrounding_loops_df["var"] == v]["boundary"].values[0]
            return result

    def _accesses_though_all_dimensions(self, reference):
        # check dimension of the array, multiply accesses based on type of expression
        surrounding_loops_df = reference["surrounding_loops"]
        count_current_refs = 1
        # for argument expressions on every dimensions
        for arg_exp in reference["arg_expressions"]:
            arg_exp = str(simplify(arg_exp)).replace(" ", "")

            if not arg_exp.isdigit():
                # TODO add case of two variables
                count_current_refs *= self._accesses_through_one_dimension(arg_exp, surrounding_loops_df)
        return count_current_refs

    @staticmethod
    def _min_expression(var_boundary, df_of_expressions):
        df_of_expressions.sort_values("value")
        if len(df_of_expressions) == 1:
            return var_boundary
        array_of_int_values = df_of_expressions["value"].values
        variations_arrays = list(variations([0, 1], len(array_of_int_values) - 1, repetition=True))
        str_formula = "min("
        for el in variations_arrays:

            indices = [index for index, element in enumerate(el) if element == 1]

            distinct = str(el.count(0) + 1) + "*" + str(var_boundary)
            overlapped = 0
            ind = None
            merging = False
            for j in range(0, len(indices)):
                if not merging:
                    ind = indices[j]
                if len(indices) > j + 1 and indices[j + 1] == indices[j] + 1:
                    merging = True
                else:
                    overlapped += abs(array_of_int_values[indices[j] + 1] - array_of_int_values[ind])
                    merging = False

            str_formula += distinct + "+" + str(overlapped)
            if variations_arrays.index(el) < len(variations_arrays) - 1:
                str_formula += ", "
            else:
                str_formula += ") "

        if str(var_boundary).isdigit():
            return eval(str_formula)
        str_formula = simplify(str_formula)
        return str_formula

    def _different_cases_of_multiplication(self, var_boundary, list_of_expressions):
        # ex. ref[2*k+2], ref[2*k+4] same since 2,4%2==0
        #     ref[2*k+3], ref[2*k+5] same since 3,5%2==1
        df_of_expressions = DataFrame(list_of_expressions)
        df_of_expressions["modulo"] = df_of_expressions.apply(lambda row: row["value"] % row["value_at_var"], axis=1)

        group_by_modulo = df_of_expressions.groupby("modulo")
        dimension_accesses_formula = 0
        for name, group in group_by_modulo:
            # value_at_var = group["value_at_var"].values[0]
            group["value"] = group.apply(lambda row: int(row["value"]), axis=1)
            dimension_accesses_formula += self._min_expression(var_boundary, group)

        return dimension_accesses_formula

    @staticmethod
    def _check_if_reference_exists(c_r, union_of_references):
        for el in union_of_references:
            if c_r["reference_name"] == el["reference_name"] and c_r["arg_expressions"] == el["arg_expressions"] and \
                    c_r["arg_variables"] == el["arg_variables"]:
                if array_equal(c_r["surrounding_loops"][["var", "boundary"]],
                               el["surrounding_loops"][["var", "boundary"]]):
                    return True
        return False


def _trivial_bound(arg_exp: str):
    """
     return variable, int multiplied with variable, int after +/- sign, sign
    """
    # a+1
    arg_exp = sub(' ', '', arg_exp)
    if bool(fullmatch("[a-zA-Z0-9_]+\+\d+(\.\d+)?", arg_exp)):
        return arg_exp.split("+")[0], None, float(arg_exp.split("+")[1]), "+", "simple"
    # a-1
    elif bool(fullmatch("[a-zA-Z0-9_]+-\d+(\.\d+)?", arg_exp)):
        return arg_exp.split("-")[0], None, float(arg_exp.split("-")[1]), "-", "simple"
    # a*2 (+ expression of integers)
    # elif bool(re.match("[a-zA-Z0-9_]+\*\d+((\+\d+)([+-\-*-/]\d+)*)*", arg_exp)):
    elif bool(fullmatch("[a-zA-Z0-9_]+\*\d+(\.\d+)?(\+\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("+")
        arg_exp.append("0")
        var_multiplication_part = arg_exp[0]
        return var_multiplication_part.split("*")[0], float(var_multiplication_part.split("*")[1]), float(arg_exp[1]), \
               "+", "mult"
    # 2*a (+ expression of integers)
    elif bool(fullmatch("\d+(\.\d+)?\*[a-zA-Z0-9_]+(\+\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("+")
        arg_exp.append("0")
        var_multiplication_part = arg_exp[0]
        return var_multiplication_part.split("*")[1], float(var_multiplication_part.split("*")[0]), float(arg_exp[1]), \
               "+", "mult"
    # a*2 (- expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+\*\d+(\.\d+)?(-\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("-")
        arg_exp.append("0")
        var_multiplication_part = arg_exp[0]
        return var_multiplication_part.split("*")[0], float(var_multiplication_part.split("*")[1]), float(arg_exp[1]), \
               "-", "mult"
    # 2*a (- expression of integers)
    elif bool(fullmatch("\d+(\.\d+)?\*[a-zA-Z0-9_]+(-\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("-")
        arg_exp.append("0")
        var_multiplication_part = arg_exp[0]
        return var_multiplication_part.split("*")[1], float(var_multiplication_part.split("*")[0]), float(arg_exp[1]), \
               "-", "mult"
    # a/2 (+ expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+/\d+(\.\d+)?((-)\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("-")
        arg_exp.append("0")
        var_division_part = arg_exp[0]
        return var_division_part.split("/")[0], (1 / float(var_division_part.split("/")[1])), float(arg_exp[1]), "+", \
               "simple"
    # a/2 (- expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+/\d+(\.\d+)?((\+)\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("+")
        arg_exp.append("0")
        var_division_part = arg_exp[0]
        return var_division_part.split("/")[0], (1 / float(var_division_part.split("/")[1])), float(arg_exp[1]), "-", \
               "simple"
    # a%2 (+ expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+%\d+(\.\d+)?(\+\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("+")
        arg_exp.append("0")
        val_modulo = arg_exp[0]
        return None, val_modulo.split("%")[1], float(arg_exp[1]), "+", "simple"
    # a%2 (- expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+%\d+(\.\d+)*(-\d+(\.\d+)?)?", arg_exp)):
        arg_exp = arg_exp.split("-")
        arg_exp.append("0")
        val_modulo = arg_exp[0]
        if val_modulo.split("%")[1].isdigit():
            return None, float(val_modulo.split("%")[1]), float(arg_exp[1]), "-", "simple"
        else:
            return None, Symbol(val_modulo.split("%")[1]), float(arg_exp[1]), "-", "simple"
    # a
    elif bool(fullmatch("[a-zA-Z0-9_]+", arg_exp)):
        return arg_exp, None, 0, "+", "simple"
    # a +- i (+- expression of integers)
    elif bool(fullmatch("[a-zA-Z0-9_]+\+\d+(\.\d+)?", arg_exp)):
        arg_exp = arg_exp.split("+")
        return arg_exp[0], None, float(arg_exp[1]), "+", "simple"
    elif bool(fullmatch("[a-zA-Z0-9_]+-\d+(\.\d+)?", arg_exp)):
        arg_exp = arg_exp.split("-|+")
        return arg_exp[0], None, float(arg_exp[1]), "-", "simple"
    elif bool(fullmatch("[a-zA-Z0-9_]+(\+[a-zA-Z0-9_]+)+([+-]\d+(\.\d+)?)*", arg_exp)):
        arg_exp = split('\\+|-', arg_exp)
        if arg_exp.__contains__("-"):
            return [arg_exp[0], arg_exp[1]], None, float(arg_exp[2]), "-", "mult_vars"
        else:
            return [arg_exp[0], arg_exp[1]], None, float(arg_exp[2]), "-", "mult_vars"
    else:
        return None, None, None, None, None
