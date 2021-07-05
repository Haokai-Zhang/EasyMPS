#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <test_exc1.py> to calculate the first excited state.
'''

from Test.params_inout import GetParams, PrintE
from AutoMPO.get_mpo import *
from AutoMPO.gen_mpo_exc1 import *
from vMPS.class_mps import mps
from EasyPlot.plot_xy import *

if __name__ == "__main__":

    # get the model parameters
    N, J, h, D, sweeps, cvg = GetParams()

    # get the Hamiltonian mpo (the variation mpo for the ground state)
    list_mpo_H = GetListHmpo(N, GetSiteHmpo_TFI(J, h))
    # get the ground state mps
    mps0 = mps(N, D, list_mpo_H)
    E0, counter = mps0.vMPS(sweeps=sweeps, cvg=cvg, if_print=True)

    # get the variation mpo for the first excited state
    list_mpo_H_excited = GetListHmpoExc1(list_mpo_H, mps0.GetListMpsData(), E0)
    # get the first excited state mps
    mps1 = mps(N, D, list_mpo_H_excited)
    E1, counter = mps1.vMPS(sweeps=sweeps, cvg=cvg, if_print=True)

    # print the results compactly
    PrintE(E1, counter)



    # measure <\sigma_z>
    list_ave_sigma_z_0 = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sz'), range(0, N))
    print('<\sigma_z> =', list_ave_sigma_z_0)

    # measure <\sigma_x>
    list_ave_sigma_x_0 = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sx'), range(0, N))
    print('<\sigma_x> =', list_ave_sigma_x_0)

    # plot <\sigma_z> and <\sigma_x>
    PlotXYn(
        range(0, N),
        [
            list_ave_sigma_z_0,
            list_ave_sigma_x_0,
        ],
        list_legend_y=[
            LatexAve('\sigma^z'),
            LatexAve('\sigma^x'),
        ],
        name_x='x',
        name_y=LatexAve('\sigma'),
        y_lim=[-1.1, 1.1],
    )


    # measure <\sigma_z>
    list_ave_sigma_z_1 = mps1.Measure1SiteOpr(2 * GenSpinOpr('Sz'), range(0, N))
    print('<\sigma_z> =', list_ave_sigma_z_1)

    # measure <\sigma_x>
    list_ave_sigma_x_1 = mps1.Measure1SiteOpr(2 * GenSpinOpr('Sx'), range(0, N))
    print('<\sigma_x> =', list_ave_sigma_x_1)

    # plot <\sigma_z> and <\sigma_x>
    PlotXYn(
        range(0, N),
        [
            list_ave_sigma_z_1,
            list_ave_sigma_x_1,
        ],
        list_legend_y=[
            LatexAve('\sigma^z'),
            LatexAve('\sigma^x'),
        ],
        name_x='x',
        name_y=LatexAve('\sigma'),
        y_lim=[-1.1, 1.1],
    )