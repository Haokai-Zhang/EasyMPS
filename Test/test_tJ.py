#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2021-07-02 16:00
Description: EasyMPS project. Run <test_tJ.py> to test <class_fsa.py>.
'''

from AutoMPO.class_fsa import fsa
from AutoMPO.class_named_data import named_data
from AutoMPO.opr_pool import GenHardCoreBosonOpr

from vMPS.class_mps import mps

if __name__ == "__main__":

    # model parameter
    Nx = 3
    Ny = 2
    t = 3.0
    J = 1.0
    # define operator with a string as name
    adagup = named_data('adagup', GenHardCoreBosonOpr('adagup'))
    aup = named_data('aup', GenHardCoreBosonOpr('aup'))
    adagdn = named_data('adagdn', GenHardCoreBosonOpr('adagdn'))
    adn = named_data('adn', GenHardCoreBosonOpr('adn'))
    F = named_data('F', GenHardCoreBosonOpr('F'))
    Id = named_data('Id', GenHardCoreBosonOpr('Id'))
    Sz = named_data('Sz', GenHardCoreBosonOpr('Sz'))
    Sp = named_data('Sp', GenHardCoreBosonOpr('Sp'))
    Sm = named_data('Sm', GenHardCoreBosonOpr('Sm'))
    ntot = named_data('ntot', GenHardCoreBosonOpr('ntot'))
    # construct finite state automata
    N = Nx * Ny
    fsa = fsa(N)
    for x in range(0, Nx):
        for y in range(0, Ny):
            # add each term in the Hamiltonian
            # horizontal bond
            if (x < Nx - 1):
                i = x * Ny + y
                j = (x + 1) * Ny + y
                fsa.Add(- t, [adagup, aup], [i, j], [Id, F, Id])
                fsa.Add(- t, [adagdn, adn], [i, j], [Id, F, Id])
                fsa.Add(- t, [aup, adagup], [i, j], [Id, F, Id])
                fsa.Add(- t, [adn, adagdn], [i, j], [Id, F, Id])
                fsa.Add(J, [Sz, Sz], [i, j])
                fsa.Add(J / 2, [Sp, Sm], [i, j])
                fsa.Add(J / 2, [Sm, Sp], [i, j])
                fsa.Add(- J / 4, [ntot, ntot], [i, j])
            # vertical bond
            if (Ny > 2) or (y < Ny - 1):
                i = x * Ny + y
                j = x * Ny + (y + 1) % Ny
                fsa.Add(- t, [adagup, aup], [i, j], [Id, F, Id])
                fsa.Add(- t, [adagdn, adn], [i, j], [Id, F, Id])
                fsa.Add(- t, [aup, adagup], [i, j], [Id, F, Id])
                fsa.Add(- t, [adn, adagdn], [i, j], [Id, F, Id])
                fsa.Add(J, [Sz, Sz], [i, j])
                fsa.Add(J / 2, [Sp, Sm], [i, j])
                fsa.Add(J / 2, [Sm, Sp], [i, j])
                fsa.Add(- J / 4, [ntot, ntot], [i, j])

    # use the constructed fsa to generate MPO
    list_mpo = fsa.GenMPO()
    # visualize MPO represented by operator symbols
    fsa.PrintSymbolMPO()
    fsa.PrintBondDimension()

    D = 15
    mps0 = mps(N, D, list_mpo)
    E0, counter = mps0.vMPS(cvg=1e-8, if_print=True, update_sites=2)