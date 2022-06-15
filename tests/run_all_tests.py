#!/usr/bin/env python3

"""@package run_all_tests
Main script to run all unit tests.
"""

from sys import exit
from unittest import TestLoader, TextTestRunner


def run():
    suite = TestLoader().loadTestsFromNames(
        [
            # Has to be done first to establish the sqlalchemy database mapping.
            'test_db_mapping',
            'test_database',
            'test_math_utils',
            'test_sample_interval',
            # Folder descriptions.
            'test_parser',
            'test_parser_util',
            # 'test_impl_skeleton',
            'test_ivp',
            'test_kernel_template',
            'test_machine',
            'test_ode_method',
            ###
            'test_performance_model',
            # 'test_train_communication',
            # 'test_utils_kerncraft',
            # 'test_benchmark',
            # Code generator.
            'test_pmodel_generator',

            # 'test_train_impl',  # TODO
            # 'test_train_kernel', # TODO
            # 'test_impl_skeleton',  # TODO
            # Folder ranking.
            # 'test_ranking',
            #####
            'test_kerncraft',
            'test_db_removal',
        ]
    )
    result = TextTestRunner(verbosity=2).run(suite)
    exit(0 if result.wasSuccessful() else 1)


if __name__ == '__main__':
    # Run the application.
    run()
