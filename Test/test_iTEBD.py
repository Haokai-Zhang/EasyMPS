#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-11-15 10:08
Description: EasyMPS project. Run <test_iTEBD.py> to test this iTEBD code.
'''

from iTEBD.class_iTEBD import *

if __name__ == "__main__":

    # initialize an instance of class iTEBD
    iTEBD = mps_iTEBD(GetPairHam_HaldanePhase(J=1,
                                              Uzz=1.5,
                                              Bx=0.,
                                              dim_spin=3),
                      D=9,
                      init_type='sz0',
                      )
    # if Uzz > 1, init_type = 'sz0'
    # if Uzz < 1, init_type = 'afm' or 'rand'

    # imaginary time evolution by logarithmic hierarchy
    iTEBD.LogEvo(max_dt=0.1,
                 max_exponent=3,
                 cvg=1e-8,
                 max_iter=1000)

    # print Schmidt weights and entanglement spectrum
    print(iTEBD.GetSchmidt(bond_head_site=0))
    #print(iTEBD.GetSchmidt(bond_head_site=1))
    print(iTEBD.GetEntangleSpec(bond_head_site=0))
    #print(iTEBD.GetEntangleSpec(bond_head_site=1))

    # print energy averaged over bonds
    print('E = %f' % iTEBD.MeasureAveE())