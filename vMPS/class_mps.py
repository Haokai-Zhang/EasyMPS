#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <class_mps.py> defines the class of matrix product state using <class_site_tensor.py>.
'''

from vMPS.class_site_tensor import site_tensor
from EasyData.data_process import *

class mps(object):
    '''
    class of mps combined with mpo
    state tensor for a single site ----contraction over sites---> matrix product state
    '''

    # initialize
    def __init__(self, N, D, list_mpo):
        '''
        :param N: number of site.
        :param D: bond dimension.
        :param list_mpo: a list of mpo corresponding to each site.
        :return: initialized mps.
        '''
        # create mps from 2 endpoints
        self.end_l = site_tensor(0, 0, 0)
        self.end_r = site_tensor(0, 0, 0)
        self.end_l.neighbor_r = self.end_r
        self.end_r.neighbor_l = self.end_l

        # create a list of site tensor
        list_mps = ([site_tensor(1, D, list_mpo[0])]
                    + [site_tensor(D, D, list_mpo[index_mpo + 1]) for index_mpo in range(N - 2)]
                    + [site_tensor(D, 1, list_mpo[N - 1])])

        # insert these state tensor into the 2 endpoints chain
        the_site_number = 0
        for the_site_tensor in list_mps:
            # insert 'the_site_tensor' on the left of 'self.end_r'.
            self.end_r.InsertNewSite(the_site_tensor, "l")
            self.end_r.neighbor_l.num = the_site_number
            the_site_number += 1

        # normalize from right to left to satisfy right canonical condition
        # so that we can start to sweep from left
        the_site_tensor = self.end_r.neighbor_l
        while (not the_site_tensor.is_end):
            the_site_tensor.SvdBondL()
            the_site_tensor = the_site_tensor.neighbor_l


    def vMPS(self, sweeps=None, cvg=None, if_print=False, update_sites=2):
        '''
        :param sweeps: times of sweeps.
        :param cvg: convergence tolerance.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :param update_sites: {1, 2} determine the number of sites to be updated in a single time.
        :return: ground state energy and sweeps/convergence counter obtained by vMPS method.
        '''
        # choose finite times of sweeps or sweeps to a given convergence tolerance
        if (sweeps == None) and (cvg == None):
            print("----Error: please set sweeps or cvg !----")
            exit(1)
        elif (sweeps != None) and (cvg != None):
            print("----Error: please only choose one of {sweeps, cvg} !----")
            exit(1)
        elif (sweeps == None) and (cvg != None):
            # fixed convergence tolerance mode
            list_E0 = []
            sweeps_counter = 0
            while (IsCvg(list_E0, cvg, check_len=2) == False):

                if (if_print == True):
                    print('')
                    print('Sweeps %d BEGIN' % sweeps_counter)

                # record the energy after each step of update
                list_E_sweep = []

                # sweep from left to right
                list_E_sweep += self.Sweep("r", if_print=if_print, update_sites=update_sites)
                list_E0.append(list_E_sweep[-1])

                # sweep from right to left
                list_E_sweep += self.Sweep("l", if_print=if_print, update_sites=update_sites)
                list_E0.append(list_E_sweep[-1])

                if (if_print == True):
                    print('Sweeps %d END' % sweeps_counter)

                sweeps_counter += 1

            return list_E0[-1], sweeps_counter
        else:
            # fixed sweep times mode
            list_E0 = []
            for sweeps_counter in range(0, sweeps):

                # print the procedure of sweep
                if (if_print == True):
                    print('')
                    print('Sweeps %d BEGIN' % sweeps_counter)

                # record the energy after each step of update
                list_E_sweep = []

                # sweep from left to right
                list_E_sweep += self.Sweep("r", if_print=if_print, update_sites=update_sites)
                list_E0.append(list_E_sweep[-1])

                # sweep from right to left
                list_E_sweep += self.Sweep("l", if_print=if_print, update_sites=update_sites)
                list_E0.append(list_E_sweep[-1])

                if (if_print == True):
                    print('Sweeps %d END' % sweeps_counter)

            cvg_counter = list_E0[-1] - list_E0[-2]

            return list_E0[-1], cvg_counter

    def Sweep(self, towards, if_print=False, update_sites=2):
        '''
        :param towards: {"l", "r"} determine the direction of sweeps.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :param update_sites: {1, 2} determine the number of sites to be updated in a single time.
        :return: a list of the lowest energies obtained by variation sweeps.
        '''
        # this function just serves as a splitter to SweepUpdate2() and SweepUpdate1()
        if (update_sites == 2):
            return self.SweepUpdate2(towards, if_print=if_print)
        elif (update_sites == 1):
            return self.SweepUpdate1(towards, if_print=if_print)
        else:
            print("----Error from Sweep(): update_sites = 1 or 2 violated !----")
            exit(1)

    def SweepUpdate2(self, towards, if_print=False):
        '''
        :param towards: {"l", "r"} determine the direction of sweeps.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :return: a list of the 2-site-updated lowest energies obtained by variation sweeps.
        '''
        list_E_sweep = []
        if (towards == "l"):
            #                           <----
            #   S---S---S---S---S---S---
            #   |   |   |   |   |   |   |   |
            #   H---H---H---H---H---H---H---H
            #   |   |   |   |   |   |   |   |
            #   S---S---S---S---S---S---
            the_site_tensor = self.end_r.neighbor_l
            while (not the_site_tensor.neighbor_l.is_end):
                list_E_sweep.append(
                    the_site_tensor.Update2Sites(towards="l",
                                                 update_with="l",
                                                 if_print=if_print,
                                                 )
                )
                the_site_tensor = the_site_tensor.neighbor_l
        elif (towards == "r"):
            #   ---->
            #        ---S---S---S---S---S---S
            #   |   |   |   |   |   |   |   |
            #   H---H---H---H---H---H---H---H
            #   |   |   |   |   |   |   |   |
            #        ---S---S---S---S---S---S
            the_site_tensor = self.end_l.neighbor_r
            while (not the_site_tensor.neighbor_r.is_end):
                list_E_sweep.append(
                    the_site_tensor.Update2Sites(towards="r",
                                                 update_with="r",
                                                 if_print=if_print,
                                                 )
                )
                the_site_tensor = the_site_tensor.neighbor_r
        else:
            print("----Error from SweepUpdate2(): towards = l or r violated !----")
            exit(1)

        return list_E_sweep

    def SweepUpdate1(self, towards, if_print=False):
        '''
        :param towards: {"l", "r"} determine the direction of sweeps.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :return: a list of the 1-site-updated lowest energies obtained by variation sweeps.
        '''
        list_E_sweep = []
        if (towards == "l"):
            #                           <----
            #   S---S---S---S---S---S---S---
            #   |   |   |   |   |   |   |   |
            #   H---H---H---H---H---H---H---H
            #   |   |   |   |   |   |   |   |
            #   S---S---S---S---S---S---S---
            the_site_tensor = self.end_r.neighbor_l
            while (not the_site_tensor.is_end):
                list_E_sweep.append(the_site_tensor.Update1Site(towards="l", if_print=if_print))
                the_site_tensor = the_site_tensor.neighbor_l
        elif (towards == "r"):
            #   ---->
            #    ---S---S---S---S---S---S---S
            #   |   |   |   |   |   |   |   |
            #   H---H---H---H---H---H---H---H
            #   |   |   |   |   |   |   |   |
            #    ---S---S---S---S---S---S---S
            the_site_tensor = self.end_l.neighbor_r
            while (not the_site_tensor.is_end):
                list_E_sweep.append(the_site_tensor.Update1Site(towards="r", if_print=if_print))
                the_site_tensor = the_site_tensor.neighbor_r
        else:
            print("----Error from SweepUpdate1(): towards = l or r violated !----")
            exit(1)

        return list_E_sweep

    def GetNorm(self):
        '''
        :return: the norm of self.
        '''
        norm_result = self.end_l.neighbor_r.GetNormR()
        # norm_result = self.end_r.neighbor_l.GetNormL()
        return norm_result.reshape(1)[0]

    def MeasureCorr(self, list_opr, list_site):
        '''
        :param list_opr: list of local operators to be measured.
        :param list_site: list of sites corresponding to these local operators.
        :return: the expectation value <ABC...> for the given mps.
        '''
        # example: list_opr = [A, B]  list_site = [2, 4], return:
        #   S---S---S---S---S---S---S---S
        #   |   |   |   |   |   |   |   |
        #   |   |   A   |   B   |   |   |
        #   |   |   |   |   |   |   |   |
        #   S---S---S---S---S---S---S---S
        list_site, list_opr = BoundSort(list_site, list_opr)
        measure_result = self.end_l.neighbor_r.GetMeasureR(list_opr, list_site)
        # measure_result = self.end_r.neighbor_l.GetMeasureL(list_opr, list_site)
        return measure_result.reshape(1)[0]

    def Measure1SiteOpr(self, opr, list_site):
        '''
        :param opr: operator to measure
        :param list_site: list of sites to measure
        :return: the expectation of opr for each sites in list_site
        '''
        list_measure_result = []
        for i in range(0, len(list_site)):
            measure_result = self.MeasureCorr([opr,],
                                              [list_site[i],],
                                              )
            list_measure_result.append(measure_result)
        return list_measure_result

    def GetEntangleSpec(self, bond_head_site, num=12):
        '''
        :param bond_head_site: the head site of the wanted bond.
        :param num: the wanted number of Schmidt weights.
        :return: entanglement spectrum corresponding to the bond after 'bond_head_site'.
        '''
        the_site_tensor = self.end_l.neighbor_r
        for the_site_num in range(0, bond_head_site):
            the_site_tensor = the_site_tensor.neighbor_r
        return (-2) * np.log(the_site_tensor.sv[0 : num])

    def GetEntangleEntropy(self):
        '''
        :return: entanglement entropy corresponding to each bond.
        '''
        list_ee = []
        the_site_tensor = self.end_l.neighbor_r
        while (not the_site_tensor.neighbor_r.is_end):
            # S = - \sum sv^2 * log(sv*2)
            sv_square = np.square(the_site_tensor.sv)
            list_ee.append(- np.dot(sv_square, np.log(sv_square)))
            the_site_tensor = the_site_tensor.neighbor_r
        return list_ee

    def GetListMpsData(self):
        '''
        :return: list of mps elements.
        '''
        list_mps_data = []
        # extract mps data from left to right
        the_site_tensor = self.end_l.neighbor_r
        while (not the_site_tensor.is_end):
            list_mps_data += [the_site_tensor.data_protected]
            the_site_tensor = the_site_tensor.neighbor_r
        return list_mps_data