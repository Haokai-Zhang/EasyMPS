#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <test_vMPS_iTEBD.py> to test entanglement spectrum obtained from vMPS and iTEBD.
'''

from EasyOpr.get_mpo import *
from vMPS.class_mps import mps
from iTEBD.class_iTEBD import *

if __name__ == "__main__":

    # set model parameters
    N = 8
    D = 6
    J = 1.
    Uzz = 0.

    # create a list of Hamiltonian matrix product operator corresponding to each site
    list_mpo_H = GetListHmpo(N, GetSiteHmpo_HaldanePhase(J=J, Uzz=Uzz, dim_spin=3))
    # create an initialized matrix product state to variate later
    mps0 = mps(N, D, list_mpo_H)
    # get the ground state by quadratic form variation
    mps0.vMPS(sweeps=None, cvg=1e-8, if_print=True, update_sites=2)
    # print entanglement spectrum obtained by vMPS
    #print(mps0.GetEntangleSpec(bond_head_site=N // 2 - 1))


    # initialize an instance of class iTEBD
    iTEBD = mps_iTEBD(GetPairHam_HaldanePhase(J=J,
                                              Uzz=Uzz,
                                              dim_spin=3),
                      D=9,
                      init_type='rand')
    # imaginary time evolution by logarithmic hierarchy
    iTEBD.LogEvo(max_dt=0.1,
                 max_exponent=3,
                 cvg=1e-9,
                 max_iter=1000)
    # print Schmidt weights and entanglement spectrum
    #print(iTEBD.GetEntangleSpec(bond_head_site=0))
    print(N * iTEBD.MeasureAveE())

