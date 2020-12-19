#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <plot_ent_entr.py> to plot entanglement entropy.
'''

from EasyOpr.get_mpo import *
from vMPS.class_mps import mps
from EasyPlot.plot_xy import *

if __name__ == "__main__":

    # number of sites
    N = 12

    # run for each h
    list_ee = []
    list_h = [0.1, 0.7, 1.0, 1.5, 10.0]
    list_h_latex = []
    for h in list_h:
        list_h_latex.append(LatexForm('h=%.1f' % h))

        # create a list of Hamiltonian matrix product operator corresponding to each site
        list_mpo_H = GetListHmpo(N=N, mpo_single_site=GetSiteHmpo_TFI(J=1, h=h))

        # create an initialized matrix product state to variate later
        mps0 = mps(N=N, D=32, list_mpo=list_mpo_H)

        # get the ground state by quadratic form variation
        E0, counter = mps0.vMPS(cvg=1e-7, if_print=True)

        # get the entanglement entropy of every bipartition
        list_ee.append(mps0.GetEntangleEntropy())

    # plot entanglement entropy
    PlotXYn(
        range(0, N-1),
        list_ee,
        name_x=LatexForm('x'),
        name_y=LatexForm('S'),
        num_column=3,
        line_width=1,
        marker_edge_width=1.2,
        y_lim=[-0.1, 0.9],
        list_legend_y=list_h_latex,
    )
    # PlotLogXY(
    #     np.sin(np.arange(1, N//2 + 1) * np.pi / N) * N / np.pi,
    #     list_ee[:N//2],
    #     name_x=LatexForm('\\tilde{x}'),
    #     name_y=LatexForm('S'),
    #     head_cut=N//6,
    #     marker_shape='-o',
    #     marker_color='red',
    #     marker_face_color='red',
    #     line_width=1,
    #     marker_edge_width=1.2,
    # )
    #
    # # measure <\sigma_z>
    # list_ave_sigma_z = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sz'), range(0, N))
    # print('<\sigma_z> =\n', list_ave_sigma_z)
    #
    # # measure <\sigma_x>
    # list_ave_sigma_x = mps0.Measure1SiteOpr(2 * GenSpinOpr('Sx'), range(0, N))
    # print('<\sigma_x> =\n', list_ave_sigma_x)