#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <test_vMPS.py> to test vMPS algorithm.
'''

from Test.params_inout import GetParams, PrintE
from EasyOpr.get_mpo import *
from vMPS.class_mps import mps
from EasyPlot.plot_xy import *

if __name__ == "__main__":

    # get the model parameters
    N, J, h, D, sweeps, cvg = GetParams()

    # create a list of Hamiltonian matrix product operator corresponding to each site
    list_mpo_H = GetListHmpo(N, GetSiteHmpo_TFI(J, h))

    # create an initialized matrix product state to variate later
    mps0 = mps(N, D, list_mpo_H)

    # get the ground state by quadratic form variation
    E0, counter = mps0.vMPS(sweeps=sweeps, cvg=cvg, if_print=True, update_sites=2)

    # print the results compactly
    PrintE(E0, counter)

    # measure <\sigma_z>
    list_ave_sigma_z = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sz'), range(0, N))
    print('<\sigma_z> =\n', list_ave_sigma_z)

    # measure <\sigma_x>
    list_ave_sigma_x = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sx'), range(0, N))
    print('<\sigma_x> =\n', list_ave_sigma_x)

    # plot <\sigma_z> and <\sigma_x>
    PlotXYn(range(0, N),
            [list_ave_sigma_z,
             list_ave_sigma_x,
             ],
            list_legend_y=[LatexAve('\sigma^z'),
                           LatexAve('\sigma^x'),
                           ],
            name_x='x',
            name_y=LatexAve('\sigma'),
            )