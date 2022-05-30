"""@package codegen.computation_dsl.computation_dsl
Definition of the computation DSL used to described computation statements used in KernelTemplate codes.

@author: Hana Shatri, Johannes Seiferth
"""

from lark import Lark, Tree

code_grammar = r"""
    ?start: lhs ASSIGN rhs

    lhs: arg_array | NAME
    rhs: expr

    expr: OP? arg_expr (OP "%"? arg_expr)*

    brackets: LB expr RB

    LB: "("
    RB: ")"

    arg_array: NAME ("[" expr "]")+       -> arg_array
    arg_expr: brackets | arg_array | NAME | INT | FLOAT

    ASSIGN: PREFIX | SUFFIX | "="
    PREFIX: (OP)? "="
    SUFFIX: "=" (OP)?

    OP: "+" | "-" | "*" | "%"| "/"

    %import common.CNAME -> NAME
    %import common.INT
    %import common.FLOAT
    %import common.WS_INLINE
    %ignore WS_INLINE
"""


def parse_lark_grammar(lark_str: str) -> Tree:
    """Park lark grammar file.

    Parameters:
    -----------
    lark_str: str
        Lark grammar file to be parsed.

    Returns:
    --------
    Lark.Tree
        Parsed Lark tree.
    """
    parser: Lark = Lark(code_grammar, parser='lalr')
    return parser.parse(lark_str)
