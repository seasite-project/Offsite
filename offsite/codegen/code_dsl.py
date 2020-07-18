"""@package code_dsl
Definition of the code DSL used to described ImplSkeleton and KernelTemplate codes.
"""

code_grammar = """
    start: (line)+ 

    line: code_line | empty_line
    ?code_line: code_keyword NEWLINE? -> code
    ?empty_line: NEWLINE

    ?code_keyword: code_pragma | code_loop_start | code_loop_end | code_computation | code_communication | code_kernel | code_pmodel | code_swp

    // Supported keywords
    code_pragma: keyword_pragma arg_word -> pragma
    ?keyword_pragma: "%PRAGMA"

    code_loop_start: keyword_loop_start arg_word arg_expression (arg_word_or_int)* -> loop_s
    ?keyword_loop_start: "%LOOP_START" -> kw_loop_s

    code_loop_end: keyword_loop_end arg_word -> loop_e
    ?keyword_loop_end: "%LOOP_END" -> kw_loop_e

    code_computation: keyword_computation arg_word -> comp
    ?keyword_computation: "%COMP"

    code_communication: keyword_communication arg_word (arg_word | arg_array)? -> comm
    ?keyword_communication: "%COM"

    code_pmodel: keyword_pmodel -> pmodel
    ?keyword_pmodel: "%PMODEL"

    code_kernel: keyword_kernel arg_word (arg_word | arg_array)? -> kernel
    ?keyword_kernel: "%KERNEL"

    code_swp: keyword_swp arg_word arg_word (arg_word | arg_pointer) -> swap
    ?keyword_swp: "%SWAP"

    // Keyword arguments
    ?arg_word: CNAME
    ?arg_word_or_int: CNAME | INT
    arg_expression: (CNAME | INT) (OPERATOR (CNAME | INT))?
    arg_array: CNAME "[" CNAME "]"
    ?arg_pointer: CNAME ("*")+

    // Custom terminals
    OPERATOR: "-"|"+"|"*"|"/"

    // Common terminals
    %import common.INT
    %import common.CNAME
    %import common.NEWLINE

    // Ignore whitespaces
    %ignore " "
"""
