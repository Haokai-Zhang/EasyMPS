#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-16 16:40
Description: EasyMPS project. Run <plot_TFI_mh.py> to plot m-h curve of TFI model using iDMRG.
'''

from iDMRG.class_iDMRG import *
from AutoMPO.get_mpo import *
from EasyPlot.plot_xy import *

if __name__ == "__main__":

    list_h = np.arange(0.02, 2, 0.04)
    list_Sx = []
    for h in list_h:
        # initialize an instance of class iDMRG
        iDMRG = mps_iDMRG(
            GetListHmpo(
                N=3,
                mpo_single_site=GetSiteHmpo_TFI(J=1, h=h),
                if_pseudo_index=False,
                ),
            D=20,
        )
        # grow sites until convergence
        E0, N = iDMRG.IterGrow(cvg=1e-5)
        print("h=%.3f  E=%.7f  N=%d" % (h, E0, N))
        # calculate expectation of spin operator
        list_Sx.append(iDMRG.Measure1Site(GenSpinOpr('Sx')))

    PlotTwin(
        list_h,
        list_Sx,
        name_x=LatexForm('h'),
        name_y1=LatexAve('S^x'),
        name_y2=LatexForm('\chi'),
    )