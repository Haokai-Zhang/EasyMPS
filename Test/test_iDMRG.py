#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-16 16:40
Description: EasyMPS project. Run <test_iDMRG.py> to test iDMRG code.
'''

from iDMRG.class_iDMRG import *
from AutoMPO.get_mpo import *
from EasyPlot.plot_xy import *

if __name__ == "__main__":

    # initialize an instance of class iDMRG
    iDMRG = mps_iDMRG(
        GetListHmpo(
            N=3,
            mpo_single_site=GetSiteHmpo_Heisenberg(J=1, h=0.),
            if_pseudo_index=False,
        ),
        D=8,
    )

    # grow sites until convergence
    list_E0, list_N = iDMRG.IterGrow(
        cvg=1e-4,
        if_print=True,
        if_return_list=True,
    )

    # plot energy v.s. number of sites
    E_exact = 0.25 - np.log(2)
    PlotXY(
        list_N,
        list_E0,
        name_x='N',
        name_y='E',
        marker_shape='-',
        dash_y=E_exact,
    )

    # plot again on log scale
    PlotLogXLogY(
        list_N,
        np.asarray(list_E0) - E_exact,
        name_x='N',
        name_y='E',
        head_cut=len(list_N)//2,
    )

    # print expectation of spin operator
    print('\n<Sx> = %.6f' % iDMRG.Measure1Site(GenSpinOpr('Sx')))
    print('<Sy> = %.6f' % iDMRG.Measure1Site(GenSpinOpr('Sy')))
    print('<Sz> = %.6f' % iDMRG.Measure1Site(GenSpinOpr('Sz')))
