"""@package code_dsl
Definition of the code DSL used to described ImplSkeleton and KernelTemplate codes.
"""

from lark import Lark, Tree
from lark.indenter import Indenter

code_grammar = r"""
    ?start: _NL* line+
    line: _kw _NL [_INDENT line+ _DEDENT]
    
    // Supported keywords
    _kw: kw_pragma | kw_loop | kw_computation | kw_communication | kw_pmodel | kw_kernel | kw_swap
    
    kw_pragma: "%PRAGMA" arg_word -> pragma
    kw_loop: "%LOOP" arg_word arg_expression (arg_word_or_int)* -> loop
    kw_computation: "%COMP" arg_word -> comp
    kw_communication: "%COM"  arg_word (arg_word | arg_array)? -> comm
    kw_pmodel: "%PMODEL" -> pmodel
    kw_kernel: "%KERNEL" arg_word (arg_word | arg_array)? -> kernel
    kw_swap: "%SWAP" arg_word arg_word (arg_word | arg_pointer) -> swap
          
    // Keyword arguments
    ?arg_word: NAME
    ?arg_word_or_int: NAME | INT
    ?arg_expression: (NAME | INT) (OPERATOR (NAME | INT))?
    ?arg_array: NAME "[" NAME "]"
    ?arg_pointer: NAME ("*")+
    
    // Custom terminals
    OPERATOR: "-"|"+"|"*"|"/"
    _NL: /(\r?\n[\t ]*)+/
    
    %import common.CNAME -> NAME
    %import common.INT
    %import common.WS_INLINE
    
    %declare _INDENT _DEDENT
    %ignore WS_INLINE
"""


class CodeIndenter(Indenter):
    NL_type = '_NL'
    OPEN_PAREN_types = []
    CLOSE_PAREN_types = []
    INDENT_type = '_INDENT'
    DEDENT_type = '_DEDENT'
    tab_len = 8


def parse_lark_grammar(lark_str: str) -> Tree:
    parser = Lark(code_grammar, parser='lalr', postlex=CodeIndenter())
    return parser.parse(lark_str)
