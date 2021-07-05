#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-11-15 10:08
Description: EasyMPS project. Run <plot_ent_spec.py> to plot entanglement spectrum vs. model parameter.
'''

from iTEBD.class_iTEBD import *
from EasyPlot.plot_xy import *
from DataToolBox.data_process import *

if __name__ == "__main__":

    # model parameters
    list_Uzz = [
        -0.7,
        0.2,
        1.5,
    ]
    # virtual bond dimension
    D = 12
    # the wanted number of Schmidt weights
    num_spec = 12
    # initialization type
    list_init_type = [
        'afm',
        'afm',
        'sz0',
    ]

    # create a list to record entanglement spectrum
    list_entangle_spec = []

    # run iTEBD for each set of model parameters
    for index_para in range(0, len(list_Uzz)):
        # print Uzz
        Uzz = list_Uzz[index_para]
        print('\nUzz = %f' % Uzz)
        # initialize an instance of class iTEBD
        iTEBD = mps_iTEBD(GetPairHam_HaldanePhase(J=1, Uzz=Uzz, Bx=0., dim_spin=3),
                          D=D,
                          init_type=list_init_type[index_para])
        # imaginary time evolution by logarithmic hierarchy
        iTEBD.LogEvo(max_dt=0.1,
                     max_exponent=3,
                     cvg=1e-8,
                     max_iter=1000)

        # print entanglement spectrum
        entangle_spec = iTEBD.GetEntangleSpec(num=num_spec)
        print(entangle_spec)
        list_entangle_spec.append(entangle_spec)

    # plot entanglement spectrum
    PlotXnYn(([ShiftDeg(list_Uzz[i], list_entangle_spec[i]) for i in range(0, len(list_Uzz))]),
             list_entangle_spec,
             name_x=LatexForm('U_{zz}/J'),
             name_y=LatexForm('-2log(\lambda_\\alpha)'),
             list_color=['black'],
             list_marker=['+'],
             marker_size=13,
             marker_edge_width=1.3,
             x_lim=[-1.0, 2.0],
             y_lim=[-0.2, 13],
             )