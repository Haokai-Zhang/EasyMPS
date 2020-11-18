#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <plot_excited.py> plot the energies of the ground state
and first excited state as a function of transverse field.
'''

from EasyOpr.get_mpo import *
from vMPS.class_mps import mps
from Test.test_excited import GetExcitedHmpo
from EasyPlot.plot_xy import *
import numpy as np

if __name__ == "__main__":

    # set the model parameters
    N, J, D, sweeps, cvg = 8, 1.0, 6, None, 1e-9

    list_h = []
    list_E0 = []
    list_E1 = []
    for h in np.arange(0, 2, 0.1):
        list_h.append(h)

        # get the Hamiltonian mpo (the variation mpo for the ground state)
        list_mpo_H = GetListHmpo(N, GetSiteHmpo_TFI(J, h))
        # get the ground state mps
        mps0 = mps(N, D, list_mpo_H)
        E0, counter = mps0.vMPS(sweeps=sweeps, cvg=cvg)

        # get the variation mpo for the first excited state
        list_mpo_H_excited = GetExcitedHmpo(list_mpo_H, mps0, E0)
        # get the first excited state mps
        mps1 = mps(N, D, list_mpo_H_excited)
        E1, counter = mps1.vMPS(sweeps=sweeps, cvg=cvg)
        list_E0.append(E0 - E0)
        list_E1.append(E1 - E0)

    # plot E0, E1
    PlotXYn(list_h,
            [list_E0,
             list_E1,
             ],
            list_ylabel=[LatexForm('E_0'),
                         LatexForm('E_1'),
                         ],
            name_x=LatexForm('g'),
            name_y=LatexForm('E_n - E_0'),
            )