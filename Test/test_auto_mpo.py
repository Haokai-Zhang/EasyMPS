#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2021-07-02 15:00
Description: EasyMPS project. Run <test_auto_mpo.py> to test <class_fsa.py>.
'''

from AutoMPO.class_fsa import fsa
from AutoMPO.class_named_data import named_data
from AutoMPO.opr_pool import GenSpinOpr

from vMPS.class_mps import mps

if __name__ == "__main__":
    # model parameter
    N = 5
    J = 1.0
    Jp = 0.5
    g = 0.5
    # define operator with a string as name
    Sz = named_data('Sz', GenSpinOpr('Sz'))
    Sx = named_data('Sx', GenSpinOpr('Sx'))
    # construct finite state automata
    fsa = fsa(N, identity_name='Id')
    for i in range(0, N):
        # add each term in the Hamiltonian
        if (i < N - 2):
            fsa.Add(Jp, [Sz, Sz], [i, i + 2], print_form='all')
        if (i < N - 1):
            fsa.Add(J, [Sz, Sz], [i, i + 1], print_form='all')
        fsa.Add(g, [Sx], [i], print_form='all')
    # use the constructed fsa to generate MPO
    list_mpo = fsa.GenMPO()
    # visualize MPO represented by operator symbols
    fsa.PrintSymbolMPO()


    D = 4
    mps0 = mps(N, D, list_mpo)
    E0, counter = mps0.vMPS(cvg=1e-9, if_print=True, update_sites=2)