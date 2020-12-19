#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <get_mpo.py> contains the Hamiltonian and other physical operators
in forms of matrix product operator.
'''

import numpy as np
from EasyOpr.gen_spin_opr import GenSpinOpr

def GetSiteHmpo_Heisenberg(J, h, dim_spin=2):
    '''
    :param J: Heisenberg coupling constant.
    :param h: magnetic field strength.
    :param dim_spin: Hilbert space dimension of single spin.
    :return: single site matrix product operator of Heisenberg model, which is a 4-order tensor. (5,5,2,2).
    '''
    # generate spin operator
    Sp, Sm, Sx, Sy, Sz, Id, S0 = GenSpinOpr(dim_spin=dim_spin)
    # H = \sum_<ij> J * kronecker(S_i, S_j) + \sum_i (-h) * S^z_i
    mpo_H_single_site_Heisenberg = np.array([[Id, S0, S0, S0, S0],
                                             [Sp, S0, S0, S0, S0],
                                             [Sm, S0, S0, S0, S0],
                                             [Sz, S0, S0, S0, S0],
                                             [- h * Sz, J / 2 * Sm, J / 2 * Sp, J * Sz, Id]])
    return mpo_H_single_site_Heisenberg
    # note that it forms a 4-order tensor, not a block matrix
    # exact ground state energy from Bethe ansatz = (1/4 - ln2) = -0.44314718

def GetSiteHmpo_TFI(J, h, dim_spin=2):
    '''
    :param J: Ising coupling constant.
    :param h: transverse field strength.
    :param dim_spin: Hilbert space dimension of single spin.
    :return: single site matrix product operator of the transverse-field Ising (TFI) model.
    '''
    # generate spin operator
    Sp, Sm, Sx, Sy, Sz, Id, S0 = GenSpinOpr(dim_spin=dim_spin)
    sigma_z = 2 * Sz
    sigma_x = 2 * Sx
    # H = (-J) * \sum_<i,i+1>( \sigma^z_i * \sigma^z_j + h * \sigma^x_i)
    mpo_H_single_site_TFI = np.array([[Id, S0, S0],
                                      [sigma_z, S0, S0],
                                      [- J * h * sigma_x, - J * sigma_z, Id]])
    return mpo_H_single_site_TFI

def GetSiteHmpo_HaldanePhase(J, Uzz, Bx=0, dim_spin=3):
    '''
    :param J: Heisenberg coupling constant.
    :params Uzz: spin suppression along z direction.
    :param Bx: transverse field strength.
    :param dim_spin: Hilbert space dimension of single spin.
    :return: single site matrix product operator of generalized spin-1 Heisenberg model, as a 4-order tensor.
    '''
    # generate spin operator
    Sp, Sm, Sx, Sy, Sz, Id, S0 = GenSpinOpr(dim_spin=dim_spin)
    # H = J * \sum{ S_i * S_j } + Uzz * \sum{ (S^z_i)^2 } + Bx * \sum{ S^x_i }
    Sz2 = np.dot(Sz, Sz)
    mpo_H_single_site_HaldanePhase = np.array([[Id, S0, S0, S0, S0],
                                               [Sp, S0, S0, S0, S0],
                                               [Sm, S0, S0, S0, S0],
                                               [Sz, S0, S0, S0, S0],
                                               [Bx * Sx + Uzz * Sz2,
                                                J / 2 * Sm, J / 2 * Sp, J * Sz, Id]])
    return mpo_H_single_site_HaldanePhase


def GetListHmpo(N, mpo_single_site, if_pseudo_index=True):
    '''
    :param N: number of sites.
    :param mpo_single_site: single site mpo of standard form.
    :param if_pseudo_index: determine if or not to append a pseudo index to the end site mpo.
    :return: a list of mpo obtained by copying the single site mpo.
    '''
    # (open boundary condition assumed)
    # sketch of mpo
    #      2         1            1
    #      |         |            |
    #  0---H---1     H---0    0---H
    #      |         |            |
    #      3         2            2

    # mpo of the left boundary site only contains the last/first row of the standard mpo
    # (for models with nearest neighbor coupling only)
    # shape = (1) + (virtual) + (physical)
    mpo_last_row = mpo_single_site[-1].copy()
    if (if_pseudo_index == True):
        shape_tuple_leftmost = (1,) + mpo_last_row.shape
        mpo_leftmost = mpo_last_row.copy().reshape(shape_tuple_leftmost)
    else:
        mpo_leftmost = mpo_last_row.copy()

    # mpo of the right boundary site only contains the last/first row of the standard mpo
    # (for models with nearest neighbor coupling only)
    # shape = (virtual) + (1,) + (physical)
    mpo_first_column = mpo_single_site[:, 0].copy()
    if (if_pseudo_index == True):
        shape_tuple_rightmost = (mpo_first_column.shape[0],) + (1,) + mpo_first_column.shape[1:]
        mpo_rightmost = mpo_first_column.copy().reshape(shape_tuple_rightmost)
    else:
        mpo_rightmost = mpo_first_column.copy()

    # copy the single site mpo of standard form
    list_mpo_H = [mpo_leftmost]
    for i in range(0, N - 2):
        list_mpo_H += [mpo_single_site.copy()]
    list_mpo_H += [mpo_rightmost]

    return list_mpo_H