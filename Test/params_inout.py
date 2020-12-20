#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <params_inout.py> contains the inputs and outputs.
'''

def GetParams():
    '''
    :return: parameters of the transverse-field Ising (TFI) model / Heisenberg model.
    '''
    # TFI: H = (-J) * \sum_<i,i+1>( \sigma^z_i * \sigma^z_j + h * \sigma^x_i)
    # Heisenberg: H = \sum_<ij> (-J) * kronecker(S_i, S_j) + \sum_<i> (-h) * S_i
    # open boundary condition
    N = 24                           # number of sites
    J = 1.0                         # coupling constant
    h = 0.01                          # magnetic field strength
    D = 4                           # maximum bond dimension of mps
    sweeps = None                   # times of sweeps
    cvg = 1e-5                      # convergence tolerance
                                    # 'sweeps' xor 'cvg' = None

    return N, J, h, D, sweeps, cvg

def PrintE(E0, counter):
    '''
    :param E0: ground state energy obtained by vMPS.
    :counter: sweeps/convergence counter.
    :return: None. Print ground state energy compactly.
    '''
    N, J, h, dim_bond_max, sweeps, cvg = GetParams()
    print('----------vMPS Results-----------')
    print('N = %d' % N)
    print('J = %f' % J, end='   ')
    print('h = %f' % h)
    print('D = %d' % dim_bond_max)
    if (sweeps == None):
        print('sweeps = %d' % counter, end='   ')
        print('cvg =', format(cvg, '.0e'))
    else:
        print('sweeps = %d' % sweeps, end='   ')
        print('cvg =', format(counter, '.0e'))
    print('E = %.10f' % E0)
    print('---------------------------------')