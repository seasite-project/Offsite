"""@package codegen.computation_dsl.computation_tree
Definition of class ComputationTreeVisitor.

@author: Hana Shatri
"""

from re import match
from typing import Dict, List

from lark import Visitor, tree as _tree

from offsite.codegen.computation_dsl.computation_dsl import parse_lark_grammar


class ComputationTreeVisitor(Visitor):

    def __init__(self):
        super().__init__()

        self.references_list: List[Dict] = list()
        self.expr_str: str = ""
        self.variables: List[str] = list()

    def get_references(self, computation_string: str) -> List[Dict[str, List[str]]]:
        """
            Returns the list of dictionaries of: reference name, arguments expressions, and arguments variables.
        """
        try:
            lark_tree = parse_lark_grammar(computation_string)
            self.visit(lark_tree)
            return self.references_list
        except:
            return self.references_list

    def arg_array(self, tree: _tree.Tree):
        reference_item: dict = {"reference_name": "", "arg_expressions": [], "arg_variables": []}

        for child in tree.children:
            if isinstance(child, _tree.Tree):
                self.expression(child)
                reference_item["arg_expressions"].append(self.expr_str)
                self.expr_str = ""
                reference_item["arg_variables"] += self.variables
                self.variables.clear()
            else:
                reference_item["reference_name"] = str(child)

        if reference_item not in self.references_list:
            self.references_list.append(reference_item)

    def expression(self, tree: _tree.Tree):
        for child in tree.children:
            if isinstance(child, _tree.Tree):
                self.expression(child)
            else:
                self.expr_str += child
                # Is the variable either a digit or operator?
                if not (match(r'^-?\d+(?:\.\d+)*$', child) or child in ["+", "-", "*", "/", "%", "(", ")"]):
                    if not self.variables.__contains__(child):
                        self.variables.append(str(child))
