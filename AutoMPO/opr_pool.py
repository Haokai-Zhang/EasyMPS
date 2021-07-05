#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <opr_pool.py> contains functions generating operators.
'''

import numpy as np

def GenSpinOpr(opr_name=None, dim_spin=2):
    '''
    :param opr_name: the name of the spin operator wanted
    :param dim_spin: Hilbert space dimension of single spin
    :return: spin operator matrix corresponding to S = (dim_spin - 1) / 2
    '''
    max_Sz = (dim_spin - 1) * 0.5
    SpSm_shoulder = [np.sqrt((2 * max_Sz - i) * (i + 1)) for i in range(0, dim_spin - 1)]
    Sp = np.diag(SpSm_shoulder, 1)
    Sm = np.diag(SpSm_shoulder, -1)
    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5j * (Sp - Sm)
    Sz = np.diag([(max_Sz - i) for i in range(0, dim_spin)])
    Id = np.identity(dim_spin)
    S0 = np.zeros((dim_spin, dim_spin))
    if (opr_name == 'Sz'):
        return Sz
    elif (opr_name == 'Sx'):
        return Sx
    elif (opr_name == 'Sy'):
        return Sy
    elif (opr_name == 'Sp'):
        return Sp
    elif (opr_name == 'Sm'):
        return Sm
    elif(opr_name == 'Id'):
        return Id
    elif (opr_name == 'real'):
        return Sz, Sp, Sm, Sx, Id
    else:
        return Sp, Sm, Sx, Sy, Sz, Id, S0

def GenHardCoreBosonOpr(opr_name=None):
    '''
    :param opr_name: the name of the boson operator wanted
    :return: boson operator matrix
    '''
    adagup = np.diag([1.], 2)
    adagdn = np.diag([0., 1.], 1)
    aup = adagup.transpose()
    adn = adagdn.transpose()
    nup = np.dot(adagup,aup)
    ndn = np.dot(adagdn,adn)
    ntot = nup + ndn
    # F = (-1)^ntot = (1 - 2*nup)(1 - 2*ndn)
    # note that 1 - np.zeros((2, 2)) = [[1,1],[1,1]]
    F = np.diag([-1., -1., 1.,])
    Id = np.identity(3)
    Sz = np.diag([0.5, -0.5, 0.])
    Sp = np.diag([1., 0.], 1)
    Sm = np.diag([1., 0.], -1)
    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5j * (Sp - Sm)
    if (opr_name == 'adagup'):
        return adagup
    elif (opr_name == 'adagdn'):
        return adagdn
    elif (opr_name == 'aup'):
        return aup
    elif (opr_name == 'adn'):
        return adn
    elif (opr_name == 'nup'):
        return nup
    elif (opr_name == 'ndn'):
        return ndn
    elif (opr_name == 'ntot'):
        return ntot
    elif (opr_name == 'F'):
        return F
    elif (opr_name == 'Id'):
        return Id
    elif (opr_name == 'Sz'):
        return Sz
    elif (opr_name == 'Sx'):
        return Sx
    elif (opr_name == 'Sy'):
        return Sy
    elif (opr_name == 'Sp'):
        return Sp
    elif (opr_name == 'Sm'):
        return Sm
    else:
        return adagup, adagdn, aup, adn, nup, ndn, ntot, F, Id, Sx, Sy, Sz