#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <gen_spin_opr.py> contains a function of generating spin operators.
'''

import numpy as np

def GenSpinOpr(which_opr=None, dim_spin=2):
    '''
    :param which_opr: choose the spin operators wanted
    :param dim_spin: Hilbert space dimension of single spin
    :return: spin operators corresponding to S = (dim_spin - 1) / 2
    '''
    max_Sz = (dim_spin - 1) * 0.5
    SpSm_shoulder = [np.sqrt((2 * max_Sz - i) * (i + 1)) for i in range(0, dim_spin - 1)]
    Sp = np.diag(SpSm_shoulder, 1)
    Sm = np.diag(SpSm_shoulder, -1)
    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5j * (Sp - Sm)
    Sz = np.diag([(max_Sz - i) for i in range(0, dim_spin)])
    Id = np.eye(dim_spin)
    S0 = np.zeros((dim_spin, dim_spin))
    if (which_opr == 'Sz'):
        return Sz
    elif (which_opr == 'Sx'):
        return Sx
    elif (which_opr == 'Sy'):
        return Sy
    elif (which_opr == 'Sp'):
        return Sp
    elif (which_opr == 'Sm'):
        return Sm
    elif (which_opr == 'real'):
        return Sz, Sp, Sm, Sx, Id
    else:
        return Sp, Sm, Sx, Sy, Sz, Id, S0