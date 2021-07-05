#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-16 16:40
Description: EasyMPS project. <class_iDMRG.py> defines the class of infinite density matrix renormalization group.
'''

import scipy.sparse.linalg
from DataToolBox.data_process import *

class mps_iDMRG(object):
    '''
    class of iDMRG object
    '''

    # initialize
    def __init__(self, list_mpo, D, if_memo_end_site=False):
        '''
        :param mpo: a list of mpo corresponding to (left end site, bulk site, right end site).
        :param D: bond dimension.
        :param if_memo_end_site: determine whether to memorize the end site mps.
        :return: initialized iDMRG tensor configuration.
        '''
        self.D = D
        self.N = 0
        self.mpo_end_l = list_mpo[0]
        self.mpo = list_mpo[1]
        self.mpo_end_r = list_mpo[2]
        self.d = self.mpo.shape[2]
        self.L_block = None
        self.R_block = None
        self.state_tensor = None
        self.sv = None
        self.Init2sites(if_memo_end_site=if_memo_end_site)

    # give current bond dimension a special name
    @property
    def D_L(self):
        return self.L_block.shape[0]
    @property
    def D_R(self):
        return self.R_block.shape[0]

    def Init2sites(self, if_memo_end_site=False):
        '''
        :param if_memo_end_site: determine whether to memorize the end site mps.
        :return: None. bootstrap tensor configuration of 2 sites.
        '''
        # record the number of sites
        self.N = 2

        # get the 2-site Hamiltonian from the given mpo
        H_2_sites = np.tensordot(
            self.mpo_end_l,
            self.mpo_end_r,
            axes=[0, 0],
        ).transpose([0, 2, 1, 3]).reshape(self.d ** 2, self.d ** 2)
        # axes=[0, 0] means
        #  1                 1
        #  |                 |
        #  H---0  <--->  0---H
        #  |                 |
        #  2                 2
        # transpose([0, 2, 1, 3]) means
        #  0      2
        #  |      |
        #  H------H
        #  |      |
        #  1      3
        # reshape(self.d ** 2, self.d ** 2) means
        #    /  \
        #   /    \
        #  0      1
        #  |      |
        #  H------H
        #  |      |
        #  2      3
        #   \    /
        #    \  /
        # in fact, it is simply H_1 * H_2

        # diagonalize the 2-site Hamiltonian
        E_complex, V_complex = scipy.sparse.linalg.eigsh(
            H_2_sites,
            1,
            which="SA",
        )
        # reshape the vector to obtain mps
        mps_2_sites = V_complex.reshape(self.d, self.d)
        # split the mps using svd
        U, S, VT, trunc_err = SvdTrunc(mps_2_sites, D_trunc=self.D)

        # define new state tensors from svd results
        state_tensor_0 = U.transpose([1, 0])
        state_tensor_1 = VT
        if (if_memo_end_site == True):
            self.end_l = state_tensor_0
            self.end_r = state_tensor_1

        # record the singular values
        self.sv = S

        # contract the mpo and mps to obtain the initialized L/R block
        self.L_block = np.tensordot(
            np.tensordot(
                state_tensor_0.conj(),
                self.mpo_end_l,
                axes=[1, 1],
            ),
            state_tensor_0,
            axes=[2, 1],
        )
        self.R_block = np.tensordot(
            np.tensordot(
                state_tensor_1.conj(),
                self.mpo_end_r,
                axes=[1, 1],
            ),
            state_tensor_1,
            axes=[2, 1],
        )
        # axes=[1, 1], [2, 1] means
        #  S---0         0---S
        #  |                 |
        #  1                 1
        #  :                 :
        #  1                 1
        #  |                 |
        #  H---0         0---H
        #  |                 |
        #  2                 2
        #  :                 :
        #  1                 1
        #  |                 |
        #  S---0         0---S

    def GetHeff2Sites(self):
        '''
        :return: the effective Hamiltonian of the 2 sites to be grown.
        '''

        # the sketch of the effective Hamiltonian H_eff
        #   .---              ---.
        #   |      |      |      |
        #   L------H------H------R
        #   |      |      |      |
        #   '---              ---'
        # as a 8-order tensor

        # contract (L---H---H---R) to get H_eff
        H_eff_order_8 = np.tensordot(
            np.tensordot(
                np.tensordot(
                    self.L_block,
                    self.mpo,
                    axes=[1, 0],
                ),
                self.mpo,
                axes=[2, 0],
            ),
            self.R_block,
            axes=[4, 1],
        ).transpose([0, 2, 4, 6, 1, 3, 5, 7])
        # axes=[1, 0] means (L---H)
        #  .---0             2
        #  |                 |
        #  L---1  <--->  0---H---1
        #  |                 |
        #  '---2             3
        # axes=[2, 0] means (L-H---H)
        #  .---0     3                  2
        #  |         |                  |
        #  L---------H---2  <--->  0 ---H--- 1
        #  |         |                  |
        #  '---1     4                  3
        # axes=[4, 1] means (L-H-H---R)
        #  .---0     2      5             0---.
        #  |         |      |                 |
        #  L---------H------H---4  <--->  1---R
        #  |         |      |                 |
        #  '---1     3      6             2---'
        # transpose([0, 2, 4, 6, 1, 3, 5, 7]) means
        #  .---0     2      4     6---.
        #  |         |      |         |
        #  L---------H------H---------R
        #  |         |      |         |
        #  '---1     3      5     7---'

        # reshape H_eff to a square matrix
        dim_H_eff_matrix = self.D_L * self.d * self.d * self.D_R
        H_eff_order_2 = H_eff_order_8.reshape(dim_H_eff_matrix, dim_H_eff_matrix)

        return H_eff_order_2

    def Grow2Sites(self, if_print=False):
        '''
        :param if_print: {True, False} determine whether to print the procedure.
        :return: ground state energy and the number of sites after growth.
        '''
        # record the growth
        self.N += 2

        # diagonalize H_eff to get the ground state
        # the eigen-vector is normalized by default of 'eigsh()'
        E_complex, V_complex = scipy.sparse.linalg.eigsh(
            self.GetHeff2Sites(),
            1,
            which="SA",
        )

        # calculate the ground state energy averaged over site
        # take the real part to avoid the infinitesimal imaginary part
        E_per_site = E_complex[0].real / self.N

        # reshape the vector to obtain mps
        mps_2_sites = V_complex.reshape(self.D_L * self.d, self.d * self.D_R)
        # split the mps using svd
        U, S, VT, trunc_err = SvdTrunc(mps_2_sites, D_trunc=self.D)

        # define new state tensors from svd results
        state_tensor_0 = U.reshape(self.D_L, self.d, len(S))
        state_tensor_1 = VT.reshape(len(S), self.d, self.D_R)

        # record the singular values and central mps
        self.sv = S
        self.state_tensor = np.tensordot(
            state_tensor_0,
            np.diag(self.sv),
            axes=[2, 0],
        )
        # axes=[2, 0] means
        #      1
        #      |
        #  0---S---2 <---> 0---sv---1

        # contract the mpo and mps to obtain the L/R block
        self.L_block = np.tensordot(
            self.L_block,
            np.tensordot(
                np.tensordot(
                    state_tensor_0.conj(),
                    self.mpo,
                    axes=[1, 2],
                ),
                state_tensor_0,
                axes=[4, 1],
            ),
            axes=[[0, 1, 2], [0, 2, 4]],
        )
        # axes=[1, 2], [4, 1] means
        #  0---S---2
        #      |
        #      1
        #      :
        #      2
        #      |
        #  0---H---1
        #      |
        #      3
        #      :
        #      1
        #      |
        #  0---S---2
        # axes=[[0, 1, 2], [0, 2, 4]]
        #  .---0  <--->  0---S---1
        #  |                 |
        #  L---1  <--->  2---H---3
        #  |                 |
        #  '---2  <--->  4---S---5
        self.R_block = np.tensordot(
            np.tensordot(
                np.tensordot(
                    state_tensor_1.conj(),
                    self.mpo,
                    axes=[1, 2],
                ),
                state_tensor_1,
                axes=[4, 1],
            ),
            self.R_block,
            axes=[[1, 3, 5], [0, 1, 2]],
        )
        # axes=[1, 2], [4, 1] means
        #  0---S---2
        #      |
        #      1
        #      :
        #      2
        #      |
        #  0---H---1
        #      |
        #      3
        #      :
        #      1
        #      |
        #  0---S---2
        # axes=[[1, 3, 5], [0, 1, 2]]
        #  0---S---1  <--->  0---.
        #      |                 |
        #  2---H---3  <--->  1---R
        #      |                 |
        #  4---S---5  <--->  2---'

        # print the procedure of growth
        if (if_print == True):
            print('N=%d' % self.N, end=' ')
            print('E=%.14f' % E_per_site, end=' ')
            print('D=%d' % self.D_L, end=' ')
            print('TruncErr=', format(trunc_err, '.1e'))

        return E_per_site, self.N

    def IterGrow(self, cvg, if_print=False, if_return_list=False):
        '''
        :param cvg: convergence tolerance.
        :param if_print: {True, False} determine whether to print the procedure.
        :param if_return_list: {True, False} determine whether to return the whole list.
        :return: (list of) ground state energy and site number.
        '''
        # record energy and number of sites
        list_E0 = []
        list_N = []
        # grow sites until convergence
        while (IsCvg(list_E0, cvg, check_len=3) == False):
            E_grow, N_grow = self.Grow2Sites(if_print=if_print)
            list_E0.append(E_grow)
            list_N.append(N_grow)
        # return list or the converged value
        if (if_return_list == True):
            return list_E0, list_N
        else:
            return list_E0[-1], list_N[-1]

    def GetEntangleSpec(self, num=8):
        '''
        :param num: the wanted number of Schmidt weights.
        :return: entanglement spectrum of self.
        '''
        return (-2) * np.log(self.sv[0: num])

    def Measure1Site(self, opr, err_imag=1e-8):
        '''
        :param opr: single site operator to be measured.
        :param err_imag: imaginary part error.
        :return: measurement result as a scalar.
        '''
        measu_result = np.tensordot(
            np.tensordot(
                self.state_tensor.conj(),
                self.state_tensor,
                axes=[[0, 2], [0, 2]],
            ),
            opr,
            axes=[[0, 1], [0, 1]],
        )
        # axes=[[0, 2], [0, 2]],[[0, 1], [0, 1]] means
        #  0---S---2
        #  :   |   :
        #  :   1   :
        #  :   :   :
        #  :   0   :
        #  :   |   :
        #  :   O   :
        #  :   |   :
        #  :   1   :
        #  :   :   :
        #  :   1   :
        #  :   |   :
        #  0---S---2
        # normalize
        measu_result /= np.square(self.sv).sum()
        # check the imaginary part
        if (measu_result.imag > err_imag):
            print("----Error from Measure1Site(): Large imaginary part.----")
            exit(1)
        else:
            return measu_result.real