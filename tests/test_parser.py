#!/usr/bin/env python3

"""@package test_parser
Tests for function of package parser.
"""

from unittest import TestCase

from offsite.descriptions.parser import ComponentDesc, ComponentDict, ComputationDesc, ComputationDict, \
    ConstantDesc, ConstantDict, DatastructDesc, DatastructDict, deserialize_obj


class TestParser(TestCase):
    def setUp(self):
        pass

    def test_types(self):
        # Test ComponentDict.
        cpnt_desc_a = ComponentDesc('memory-bound', 'x[j]=4', 'n', '0')
        cpnt_desc_b = ComponentDesc('compute-bound', 'y[j]=1', 1, 0)
        cpnt_dict = ComponentDict()
        cpnt_dict['a'] = cpnt_desc_a
        cpnt_dict['b'] = cpnt_desc_b
        cpnt_dict['a'] = cpnt_desc_a
        self.assertEquals(cpnt_desc_a, cpnt_dict.get_component('a'))
        with self.assertRaises(RuntimeError):
            cpnt_dict.get_component('c')
        cpnt_dict_str = cpnt_dict.serialize()
        cpnt_dict_cpy = deserialize_obj(cpnt_dict_str)
        self.assertEquals(cpnt_dict, cpnt_dict_cpy)
        # Test ComputationDict.
        comp_desc_a = ComputationDesc('x[j]=4')
        comp_desc_b = ComputationDesc('y[j]=1')
        comp_dict = ComputationDict()
        comp_dict['a'] = comp_desc_a
        comp_dict['b'] = comp_desc_b
        comp_dict['a'] = comp_desc_a
        self.assertEquals(comp_desc_a, comp_dict.get_computation('a'))
        with self.assertRaises(RuntimeError):
            comp_dict.get_computation('c')
        comp_dict_str = comp_dict.serialize()
        comp_dict_cpy = deserialize_obj(comp_dict_str)
        self.assertEquals(comp_dict, comp_dict_cpy)
        # Test ConstantDict.
        const_desc_a = ConstantDesc('double', 4.21)
        const_desc_b = ConstantDesc('int', 2)
        const_desc_c = ConstantDesc('double', 'A*2')
        const_dict = ConstantDict()
        const_dict['a'] = const_desc_a
        const_dict['b'] = const_desc_b
        const_dict['a'] = const_desc_a
        self.assertEquals(const_desc_a, const_dict.get_constant('a'))
        with self.assertRaises(RuntimeError):
            const_dict.get_constant('c')
        const_dict['c'] = const_desc_c
        const_dict_str = const_dict.serialize()
        const_dict_cpy = deserialize_obj(const_dict_str)
        self.assertEquals(const_dict, const_dict_cpy)
        # Test DatastructDict.
        data_desc_a = DatastructDesc('float', 1, 'n', False)
        data_desc_b = DatastructDesc('int', 2, ['s', 'n'], False)
        data_desc_c = DatastructDesc('int', 3, ['x', 'y', 'z'], True)
        data_desc_d = DatastructDesc('float', 0, '', False)
        data_dict = DatastructDict()
        data_dict['a'] = data_desc_a
        data_dict['b'] = data_desc_b
        data_dict['a'] = data_desc_a
        self.assertEquals(data_desc_a, data_dict.get_datastruct('a'))
        with self.assertRaises(RuntimeError):
            data_dict.get_datastruct('e')
        data_dict_str = data_dict.serialize()
        data_dict_cpy = deserialize_obj(data_dict_str)
        self.assertEquals(data_dict, data_dict_cpy)
        self.assertEquals(data_desc_a.dimensions_to_string(), '[n]')
        self.assertEquals(data_desc_b.dimensions_to_string(), '[s][n]')
        self.assertEquals(data_desc_c.dimensions_to_string(), '[x][y][z]')
        self.assertEquals(data_desc_d.dimensions_to_string(), '')
