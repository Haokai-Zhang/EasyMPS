#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <test_iDMRG_iTEBD.py> to test entanglement spectrum obtained from iDMRG and iTEBD.
'''

from EasyOpr.get_mpo import *
from iDMRG.class_iDMRG import *
from iTEBD.class_iTEBD import *

if __name__ == "__main__":

    # set model parameters
    D = 8
    J = 1.
    Uzz = 1.5

    # initialize an instance of class iDMRG
    iDMRG = mps_iDMRG(
        GetListHmpo(
            N=3,
            mpo_single_site=GetSiteHmpo_HaldanePhase(J=J, Uzz=Uzz, dim_spin=3),
            if_pseudo_index=False,
        ),
        D=D,
    )
    # grow sites until convergence
    E0, N = iDMRG.IterGrow(cvg=1e-5)
    print('\niDMRG')
    print('E = %.10f' % E0)
    print(iDMRG.GetEntangleSpec())

    # initialize an instance of class iTEBD
    iTEBD = mps_iTEBD(GetPairHam_HaldanePhase(J=J, Uzz=Uzz, dim_spin=3),
                      D=D,
                      init_type='rand')
    # imaginary time evolution by logarithmic hierarchy
    iTEBD.LogEvo(max_dt=0.1,
                 max_exponent=3,
                 cvg=1e-9,
                 max_iter=1000)
    # print Schmidt weights and entanglement spectrum
    print('\niTEBD')
    print('E = %.10f' % iTEBD.MeasureAveE())
    print(iTEBD.GetEntangleSpec(bond_head_site=0))