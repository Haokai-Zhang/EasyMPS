#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <gen_mpo_exc1.py> contains functions to get the variation mpo
 for the first excited state.
'''

from DataToolBox.data_process import *

def GetListHmpoExc1(list_mpo_H, list_mps0_data, E0, if_pseudo_index=True):
    '''
    :param list_mpo_H: mpo list of the original Hamiltonian.
    :param list_mps0_data: list of mps elements.
    :param E0: ground state energy.
    :param if_pseudo_index: {True, False} determine whether the pseudo-indices are used.
    :return: the variation mpo for the first excited state.
    '''
    # H' = H - E0 * | 0 > < 0 |

    # construct the ground state projector mpo by tensor product for each site
    # note that this projector has bond dimension of D^2
    # which brings us contraction process of D^6 (ground state only D^4)
    list_mpo_projector = []
    for mps_data in list_mps0_data:
        # tensor product
        direct_prod = np.tensordot(
            mps_data,
            mps_data.conj(),
            axes=0,
        )
        # axes=0 means tensor product
        #       1                          1
        #       |                          |
        #   0---S---2     tensor       0---S---2
        #       :        ========>         |
        #   0---S---2     product      3---S---5
        #       |                          |
        #       1                          4
        #
        #   1             1
        #   |             |
        #   S---0         S---0
        #   :      ====>  |
        #   S---0         S---2
        #   |             |
        #   1             3

        # combine the virtual dimensions of the projector
        # to make it consistent with the Hamiltonian mpo
        if (if_pseudo_index == False) and (mps_data.ndim == 2):
            mpo_projector = CombineIndex(
                direct_prod,
                list_index_to_combine=[[0, 2], [1, ], [3, ]],
            )
            # list_combine_pair=[[0, 2], [1, ], [3, ]] means
            #   1             1
            #   |             |
            #   S---0         S---.             1
            #   |             |    \            |
            #   |      ====>  |     --0   ===>  P---0
            #   |             |    /            |
            #   S---2         S---'             2
            #   |             |
            #   3             3
        else:
            mpo_projector = CombineIndex(
                direct_prod,
                list_index_to_combine=[[0, 3], [2, 5], [1, ], [4, ]],
            )
            # list_combine_pair=[[0, 3], [2, 5], [1, ], [4, ]] means
            #       1                         2
            #       |                         |
            #   0---S---2                 .---S---.                 2
            #       |                    /    |    \                |
            #       |       =====>    0--     |     --1   ===>  0---P---1
            #       |                    \    |    /                |
            #   3---S---5                 '---S---'                 3
            #       |                         |
            #       4                         3
        list_mpo_projector.append(mpo_projector)


    # projector multiplied by scalar (-E0)
    list_mpo_projector[0] = (- E0) * list_mpo_projector[0]

    # construct the variation mpo for the first excited state by
    # adding the mpo of Hamiltonian and projector on the virtual bond,
    # H' = H - E0 * | 0 > < 0 |
    # where the additivity is guaranteed by the same physical bond dimension
    list_mpo_H_excited = []
    for site_index in range(0, len(list_mpo_H)):
        list_mpo_H_excited.append(
            ExtendAdd(
                list_mpo_H[site_index],
                list_mpo_projector[site_index],
            )
        )
        print(list_mpo_H_excited[-1].shape)

    return list_mpo_H_excited