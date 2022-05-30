"""@package util.working_sets.count_references
Definitions of class CountReferences.

@author: Hana Shatri
"""

from re import fullmatch, match, split
from typing import List

from numpy import array_equal
from pandas import DataFrame, concat
from sympy import Symbol, simplify
from sympy.parsing.sympy_parser import parse_expr

from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType
from offsite.descriptions.parser import DatastructDict
from offsite.util.working_sets.direct_formula_definition import DirectFormulaDefinition


class CountReferences:
    def __init__(self):
        # the simplified dataframe object of computation and loop nodes that will be used through the methods
        self.code_tree_in_df = DataFrame(columns=["code_row", "indent", "type", "references", "var", "start",
                                                  "boundary", "surrounding_loops"])
        self.surrounding_loops_of_comps = DataFrame(columns=["code_row", "indent", "var", "boundary"])
        self.code_row = 0
        self.working_sets = list()

    def compute_working_sets(self, code_root: CodeNode, codegen) -> List[str]:
        self.working_sets = list()

        # return CodeTree, with all the information for loops and computation strings
        self.visit_top_down(code_root, codegen)

        self.recursive_count_of_references(self.code_tree_in_df)

        return self.working_sets

    def visit_top_down(self, node: CodeNode, codegen, end_at_nxt_pmodel_node=False):
        """ Visit the returned tree bottom up and populate the code_tree_in_df.
        """
        if node.type == CodeNodeType.PMODEL:
            if end_at_nxt_pmodel_node:
                return
            else:
                end_at_nxt_pmodel_node = True

        def _loop2df(row, indent, variable, boundary) -> DataFrame:
            return DataFrame([{"code_row": row, 'indent': indent, "var": variable, "boundary": boundary}])

        def _comp2df(row, indent, type_, references, variable, start, boundary) -> DataFrame:
            return DataFrame([{"code_row": row, "indent": indent, "type": type_, "references": references,
                               "var": variable, "start": start, "boundary": boundary}])

        #
        if node.child:
            if not self.surrounding_loops_of_comps.empty:
                indexNames = self.surrounding_loops_of_comps[
                    (self.surrounding_loops_of_comps['indent'] >= node.child.indent)].index
                self.surrounding_loops_of_comps.drop(indexNames, inplace=True)

            # append object
            if node.child.type == CodeNodeType.COMPUTATION:
                if node.child.isIVPdependent:
                    if 'RHS_input' in codegen and 'RHS_butcher_nodes' in codegen:
                        for vname, var in DatastructDict.from_data(
                                ['double {}[j]'.format(codegen['RHS_input']),
                                 'double {}'.format(codegen['RHS_butcher_nodes'])]).items():
                            node.child.references.append(
                                {'reference_name': vname, 'arg_expressions': var.size, 'arg_variables': var.size})
                    else:
                        assert False

                # each reference should have its surrounding loops, since they are treated in atomic
                for reference in node.child.references:
                    surrounding_loops_of_this_reference = self.surrounding_loops_of_comps[
                        self.surrounding_loops_of_comps['var'].isin(reference["arg_variables"])]
                    reference["surrounding_loops"] = surrounding_loops_of_this_reference

                self.code_tree_in_df = concat([self.code_tree_in_df,
                                               _comp2df(self.code_row, node.child.indent, node.child.type,
                                                        node.child.references, None, None, None)], ignore_index=True)

            elif node.child.type == CodeNodeType.LOOP:
                self.code_tree_in_df = concat([self.code_tree_in_df,
                                               _comp2df(self.code_row, node.child.indent, node.child.type, None,
                                                        node.child.var, node.child.start, node.child.boundary)],
                                              ignore_index=True)

                boundary_calc = Symbol(node.child.boundary) - int(node.child.start)
                self.surrounding_loops_of_comps = concat([self.surrounding_loops_of_comps,
                                                          _loop2df(self.code_row, node.child.indent, node.child.var,
                                                                   boundary_calc)], ignore_index=True)

            if node.child.type is not CodeNodeType.PMODEL:
                self.code_row += 1
            self.visit_top_down(node.child, codegen, end_at_nxt_pmodel_node)

        if node.next:
            if not self.surrounding_loops_of_comps.empty:
                indexNames = self.surrounding_loops_of_comps[
                    (self.surrounding_loops_of_comps['indent'] >= node.next.indent)].index
                self.surrounding_loops_of_comps.drop(indexNames, inplace=True)

            # append object
            if node.next.type == CodeNodeType.COMPUTATION:
                if node.next.isIVPdependent:
                    if 'RHS_input' in codegen and 'RHS_butcher_nodes' in codegen:
                        for vname, var in DatastructDict.from_data(
                                ['double {}[j]'.format(codegen['RHS_input']),
                                 'double {}'.format(codegen['RHS_butcher_nodes'])]).items():
                            node.next.references.append(
                                {'reference_name': vname, 'arg_expressions': var.size, 'arg_variables': var.size})
                    else:
                        assert False

                for reference in node.next.references:
                    surrounding_loops_of_this_reference = self.surrounding_loops_of_comps[
                        self.surrounding_loops_of_comps['var'].isin(reference["arg_variables"])]
                    reference["surrounding_loops"] = surrounding_loops_of_this_reference

                self.code_tree_in_df = concat([self.code_tree_in_df,
                                               _comp2df(self.code_row, node.next.indent, node.next.type,
                                                        node.next.references, None, None, None)], ignore_index=True)

            elif node.next.type == CodeNodeType.LOOP:
                self.code_tree_in_df = concat([self.code_tree_in_df,
                                               _comp2df(self.code_row, node.next.indent, node.next.type, None,
                                                        node.next.var, node.next.start, node.next.boundary)],
                                              ignore_index=True)

                boundary_calc = Symbol(node.next.boundary) - int(node.next.start)
                self.surrounding_loops_of_comps = concat([self.surrounding_loops_of_comps,
                                                          _loop2df(self.code_row, node.next.indent, node.next.var,
                                                                   boundary_calc)], ignore_index=True)

                self.code_row += 1
                self.visit_top_down(node.next, codegen, end_at_nxt_pmodel_node)

    def recursive_count_of_references(self, code_tree_df: DataFrame):
        """ Find working sets recursively for all loop nests
        while checking if there are references that refer to a non-existent outer loop
        """
        min_indent = min(code_tree_df["indent"])
        max_indent = max(code_tree_df["indent"])

        if max_indent - min_indent > 0:
            # loops for which the working set will be found
            row_values: list = list(
                code_tree_df[(code_tree_df["indent"] == min_indent) & (code_tree_df["type"] == CodeNodeType.LOOP)][
                    "code_row"].values)
            row_values.append(code_tree_df.iloc[-1]['code_row'] + 1)

            for index in range(0, len(row_values) - 1):
                chunked_code_tree = code_tree_df[(code_tree_df["code_row"] >= row_values[index]) &
                                                 (code_tree_df["code_row"] < row_values[index + 1])].copy()

                chunked_code_tree = self.handle_ref_dimensions_that_use_nonexistent_loops(chunked_code_tree)
                self.count_references(chunked_code_tree)
                # remove the outer loop, outer comp, and comps at one level deeper
                indexNames = chunked_code_tree[(chunked_code_tree['indent'] == min_indent) |
                                               ((chunked_code_tree["type"] == CodeNodeType.COMPUTATION) &
                                                (chunked_code_tree["indent"] <= min_indent + 1))].index
                chunked_code_tree.drop(indexNames, inplace=True)
                if not chunked_code_tree.empty:
                    self.recursive_count_of_references(chunked_code_tree)

    def count_references(self, code_tree_df: DataFrame):
        """ Count references by using the strategy that fits with the structure of the argument expressions

        Direct creation of formulas for:
            unique references
            non-unique references with simple bounds

        Simulation:
            for non-unique references with difficult bounds and for small n
        Approximation:
              other
        """
        unique, trivial, unique_or_simple_references, difficult_references = \
            self.check_reference_uniqueness_and_expressions_format(code_tree_df)

        dfd = DirectFormulaDefinition(code_tree_in_df=code_tree_df)

        if unique:
            formula = dfd.unique_references_count()
        elif trivial:
            formula = dfd.non_unique_simple_references_count()
        else:
            assert False

        self.working_sets.append(str(formula))

    def check_reference_uniqueness_and_expressions_format(self, code_tree_df: DataFrame):
        unique_or_simple_references = []
        difficult_references = []
        # union all references
        union_of_distinct_references = []

        for index, row in code_tree_df.iterrows():

            if row["type"] == CodeNodeType.COMPUTATION:
                comp_references = row["references"]
                # for each reference of the computation node, ex. comp node -> (A[i] = B[i][k]) where A,B are references
                for reference in comp_references:
                    # if two references are all identical, only one of them is needed for calculation
                    if not self.check_if_reference_exists(reference, union_of_distinct_references):
                        union_of_distinct_references.append(reference)

        # check if there are two or more references with the same name, (uniqueness)
        # meaning that their access expressions are different because the identical ones
        # are filtered out from previous code block

        groups = DataFrame(union_of_distinct_references).groupby("reference_name")
        if len(groups) == len(union_of_distinct_references):
            # all references are unique
            # unique: true, format -> not important
            return True, None, None, None
        else:
            are_all_references_simple = True
            # references are not all unique, so check the format of arguments expressions of non-unique references,
            # to know if the trivial bounds can be used or not
            # this is done using regex
            # * formulas of the difficult cases can then be found by simulation or Lagrange interpolation

            # iterate over groups of references with the same name
            # for each group of reference if they are of the form '...' trivial bounds can be used
            for name, group in groups:

                unique_or_simple_references.append(group["reference_name"].values[0])
                boundary_and_int_values_of_first = []
                # skip unique references
                if len(group) > 1:
                    group_element_index = 0
                    # for each reference with the same name
                    for index, reference in group.iterrows():
                        surrounding_loops_df = reference["surrounding_loops"]
                        if not surrounding_loops_df.empty:

                            reference_expression_index = 0
                            # for expression in each dimension of the reference,
                            # ex: A[expression in dimension 1][expression in dimension 2]..[expression in dimension n]
                            for exp in reference["arg_expressions"]:

                                exp = str(simplify(exp)).replace(" ", "")
                                # checking simplicity of non-unique references
                                simple, var, int_value = _check_simplicity(exp)

                                # if at least one of the expressions in the reference (e.x. ref[exp][exp]..[exp])
                                # is not simple, then the whole reference group is treated as difficult
                                if not simple:
                                    # unique: false, format: not simple
                                    are_all_references_simple = False
                                    if reference["reference_name"] in unique_or_simple_references:
                                        unique_or_simple_references.remove(reference["reference_name"])
                                    if reference["reference_name"] not in difficult_references:
                                        difficult_references.append(reference["reference_name"])

                                if are_all_references_simple:

                                    # here we check for the non-unique cases if the value attached to the variable is
                                    # the same for all. Ex. ref[2*i+1] and ref[2*i-10] are okay, but not ref[2*i], ref[3*i]
                                    # and if the variable goes up to the same bound
                                    if group_element_index == 0:
                                        # Add the var value and the int attached to it
                                        variable = \
                                            surrounding_loops_df[surrounding_loops_df["var"] == var][
                                                "boundary"].values[
                                                0]
                                        boundary_and_int_values_of_first.append([variable, int_value])
                                    else:
                                        variable = \
                                            surrounding_loops_df[surrounding_loops_df["var"] == var][
                                                "boundary"].values[
                                                0]
                                        if [variable, int_value] != boundary_and_int_values_of_first[
                                            reference_expression_index]:
                                            are_all_references_simple = False
                                            unique_or_simple_references.remove(reference["reference_name"])
                                            if reference["reference_name"] not in difficult_references:
                                                difficult_references.append(reference["reference_name"])
                                            break
                                reference_expression_index += 1
                            group_element_index += 1

            return False, are_all_references_simple, unique_or_simple_references, difficult_references

    @staticmethod
    def check_if_reference_exists(c_r, union_of_references) -> bool:
        for el in union_of_references:
            if c_r["reference_name"] == el["reference_name"] and c_r["arg_expressions"] == el["arg_expressions"] and \
                    c_r["arg_variables"] == el["arg_variables"]:
                if array_equal(c_r["surrounding_loops"][["var", "boundary"]],
                               el["surrounding_loops"][["var", "boundary"]]):
                    return True
        return False

    @staticmethod
    def handle_ref_dimensions_that_use_nonexistent_loops(tree_nodes: DataFrame) -> DataFrame:
        comp_nodes_df: DataFrame = tree_nodes[tree_nodes['type'] == CodeNodeType.COMPUTATION]
        loop_nodes_df: DataFrame = tree_nodes[tree_nodes['type'] == CodeNodeType.LOOP]

        # for every computation node
        for index, df_row in comp_nodes_df.iterrows():

            # for every reference in that computation node
            for ref in df_row["references"]:
                # remove surrounding loops that doesn't exits
                drop_indexes = []

                # iterate over surrounding loops of this reference
                for ind, s_loop in ref["surrounding_loops"].iterrows():
                    if not s_loop["code_row"] in loop_nodes_df["code_row"]:
                        drop_indexes.append(ind)
                ref["surrounding_loops"] = ref["surrounding_loops"].drop(drop_indexes)  # = temp
                # remove argument variables referring to non_existent loops --> intersection
                intersection = list(set(ref["surrounding_loops"]["var"].values) & set(ref["arg_variables"]))
                # convert argument expressions referring to that variable to one
                # simplify
                difference = list(set(ref["arg_variables"]) - set(intersection))
                difference = [Symbol(x) for x in difference]
                dict_exp_to_val = {}
                for el in difference:
                    dict_exp_to_val[el] = 1

                temp_list_of_arg_exp = ref["arg_expressions"]
                for idx, arg_expression in enumerate(ref["arg_expressions"]):
                    temp_list_of_arg_exp[idx] = str(parse_expr(arg_expression, evaluate=True).subs(dict_exp_to_val))
                ref["arg_expressions"] = temp_list_of_arg_exp
                ref["arg_variables"] = intersection
        return tree_nodes


def _check_simplicity(expression):
    # check if it is one of the simple formats
    expression = expression.replace(" ", "")
    if fullmatch("[a-zA-Z0-9_]+[*%/]\d+(\.\d+)?([+-]\d+(\.\d+)?)+", expression):
        exp_array = split("\\*|/|%|\\+|-", expression)
        var = exp_array[0]
        int_value = exp_array[1]
        return True, var, int(int_value)
    elif match("\d+(.\d+)?[*][a-zA-Z0-9_]+([+-]\d+(.\d+)?)+", expression):
        exp_array = split("\\*|\\+|-", expression)
        var = exp_array[1]
        int_value = exp_array[0]
        return True, var, int(int_value)
    elif match("[a-zA-Z0-9_]+([+-]\d+(.\d+)?)*", expression):
        exp_array = split("\\+|-", expression)
        var = exp_array[0]
        return True, var, None

    return False, None, None
