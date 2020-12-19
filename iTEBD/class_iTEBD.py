#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-11-15 10:08
Description: EasyMPS project. <class_iTEBD.py> defines a class of mps used in the iTEBD algorithm.
'''

import numpy as np
from EasyOpr.gen_spin_opr import GenSpinOpr


class mps_iTEBD(object):
    '''
    class of mps used in the iTEBD algorithm.
    '''

    # state tensor for a single site: 3-order tensor (2 virtual + 1 physical). Marked as "S" for "State".
    #             shape[1]                  (d)
    #               |                        |
    #   shape[0] ---S--- shape[2]      (D_max)---S---(D_max)

    # iTEBD configuration
    #           |          |
    #         expH--------expH
    #           |          |
    #  ---v10---S0---v01---S1---v10---
    # or
    #           |          |
    #         expH--------expH
    #           |          |
    #  ---v01---S1---v10---S0---v01---

    # initialize
    def __init__(self, H_pair, D, init_type='rand'):
        '''
        :param H_pair: two-site interaction unit of Hamiltonian.
        :param D: maximum virtual bond dimension of self.
        :param init_type: type of initial state.
        :return: initialized tensors of iTEBD configuration.
        '''

        # bond dimension
        self.d_phys = H_pair.shape[0]
        self.D_max = D
        if (init_type == 'rand'):
            # random initialization of state tensor
            self.S0 = np.random.rand(D, self.d_phys, D)
            self.S1 = np.random.rand(D, self.d_phys, D)
            # random initialization of singular values
            self.v01 = np.random.rand(D)
            self.v10 = np.random.rand(D)
        elif (init_type == 'sz0'):
            # Sz=0 initialization of state tensor
            half_sz = self.d_phys // 2
            self.S0 = np.zeros((1, self.d_phys, 1))
            self.S1 = np.zeros((1, self.d_phys, 1))
            self.S0[0, half_sz, 0] = 1.0
            self.S1[0, half_sz, 0] = 1.0
            # direct product state initialization of singular values
            self.v01 = np.ones(1)
            self.v10 = np.ones(1)
        elif (init_type == 'afm'):
            # anti-ferromagnetic initialization of state tensor
            self.S0 = np.zeros((1, self.d_phys, 1))
            self.S1 = np.zeros((1, self.d_phys, 1))
            self.S0[0, 0, 0] = 1.0
            self.S1[0, self.d_phys - 1, 0] = 1.0
            # direct product state initialization of singular values
            self.v01 = np.ones(1)
            self.v10 = np.ones(1)
        else:
            print("----Error from mps_iTEBD(): init_type = {rand, sz0, afm} !")

        # record Hamiltonian
        self.H_pair = H_pair
        self.exp_Ht = None

    # give bond dimension a special name
    @property
    def D_01(self):
        return self.S0.shape[2]

    @property
    def D_10(self):
        return self.S0.shape[0]

    def ContractState(self, head_site):
        '''
        :param head_site: the head site for contraction.
        :return: contracted state tensor.
        '''
        if (head_site == 0):
            #                  1                        1
            #                  |                        |
            # 0---v10---1  0---S0---2  0---v01---1  0---S1---2  0---v10---1
            state_tensor = np.tensordot(
                np.tensordot(
                    np.tensordot(np.diag(self.v10),
                                 self.S0,
                                 axes=[1, 0],
                                 ),
                    np.diag(self.v01),
                    axes=[2, 0],
                ),
                np.tensordot(self.S1,
                             np.diag(self.v10),
                             axes=[2, 0],
                             ),
                axes=[2, 0],
            )
        else:
            #                  1                        1
            #                  |                        |
            # 0---v01---1  0---S1---2  0---v10---1  0---S0---2  0---v01---1
            state_tensor = np.tensordot(
                np.tensordot(
                    np.tensordot(np.diag(self.v01),
                                 self.S1,
                                 axes=[1, 0],
                                 ),
                    np.diag(self.v10),
                    axes=[2, 0],
                ),
                np.tensordot(self.S0,
                             np.diag(self.v01),
                             axes=[2, 0],
                             ),
                axes=[2, 0],
            )
        return state_tensor

    def SingleEvo(self, bond_head_site):
        '''
        :param bond_head_site: the head site of the evolved bond.
        :return: None. Evolve a single bond.
        '''
        # contract the mps to be evolved
        mps_old = self.ContractState(head_site=bond_head_site)
        # rename the outside bond
        if (bond_head_site == 0):
            D_out = self.D_10
        else:
            D_out = self.D_01

        # implement exp(-Hdt) to mps
        mps_new = np.tensordot(mps_old,
                               self.exp_Ht,
                               axes=[[1, 2], [2, 3]],
                               ).transpose([0, 2, 3, 1])
        # axes=[[1, 2], [2, 3]] means
        #           0         1
        #           |         |
        #         expH-------expH
        #           |         |
        #           2         3
        #           :         :
        #           1         2
        #           |         |
        #   0---v---S----v----S---v---3
        # transpose([0, 2, 3, 1]) means
        #           2         3
        #           |         |
        #   0---v---S----v----S---v---1

        # reshape to a square matrix and truncate by svd
        U, S, VT = np.linalg.svd(mps_new.reshape(D_out * self.d_phys,
                                                 self.d_phys * D_out),
                                 full_matrices=False,
                                 )

        # the bond dimension between U&S*VT / U*S&VT might be larger than D_max due to D_max*d>D_max
        # truncate
        if (len(S) > self.D_max):
            U = U[:, 0: self.D_max]
            S = S[0: self.D_max]
            VT = VT[0: self.D_max, :]
        D_new = len(S)

        # normalization
        S /= np.sqrt(np.sum(S ** 2))

        # update self
        # print(mps_new.shape)
        if (bond_head_site == 0):
            self.v01 = S
            self.S0 = np.tensordot(
                np.diag(1 / self.v10),
                U.reshape(self.D_10, self.d_phys, D_new),
                axes=[1, 0],
            )
            self.S1 = np.tensordot(
                VT.reshape(D_new, self.d_phys, self.D_10),
                np.diag(1 / self.v10),
                axes=[2, 0],
            )
        else:
            self.v10 = S
            self.S1 = np.tensordot(
                np.diag(1 / self.v01),
                U.reshape(self.D_01, self.d_phys, D_new),
                axes=[1, 0],
            )
            self.S0 = np.tensordot(
                VT.reshape(D_new, self.d_phys, self.D_01),
                np.diag(1 / self.v01),
                axes=[2, 0],
            )

    def IterEvo(self, cvg, max_iter, if_print=False):
        '''
        :param cvg: convergence tolerance.
        :param max_iter: maximum times of iteration.
        :param if_print: determine whether to print the procedure.
        :return: None. Iterative evolution.
        '''
        # take the maximum singular value as a criterion of convergence
        max_v01 = 1
        max_v10 = 1
        # start iteration
        cvg_err = None
        for index_iter in range(0, max_iter):
            # single evolution
            self.SingleEvo(bond_head_site=0)
            self.SingleEvo(bond_head_site=1)
            # convergence error = the change of maximum singular value
            cvg_err = abs(self.v01[0] - max_v01) + abs(self.v10[0] - max_v10)
            max_v01 = self.v01[0]
            max_v10 = self.v10[0]
            # if converge, break
            if (cvg_err < cvg):
                if (if_print == True):
                    print('cvg =', format(cvg_err, '.0e'), end='    ')
                    print('iter = %d' % index_iter)
                return
        # print convergence tolerance
        if (if_print == True):
            print('cvg =', format(cvg_err, '.0e'), end='    ')
            print('iter = %d' % max_iter)

    def LogEvo(self, max_dt, max_exponent, cvg, max_iter, if_print=False):
        '''
        :param max_dt: maximum time evolution step.
        :param max_exponent: maximum exponent of the time fineness.
        :param cvg: convergence tolerance.
        :param max_iter: maximum times of iteration.
        :param if_print: determine whether to print the procedure.
        :return: None. Imaginary time evolution by logarithmic hierarchy.
        '''
        # for each hierarchy
        for exponent in range(0, max_exponent):
            # e.g. dt = 0.1, 0.01, 0.001, ...
            dt = max_dt * (0.1 ** exponent)
            if (if_print == True):
                print('dt = %f' % dt, end='    ')
            self.exp_Ht = GetExpHamPair(self.H_pair, dt)
            # iterative evolution for each dt
            self.IterEvo(cvg, max_iter, if_print)

    def GetEntangleSpec(self, num=8, bond_head_site=0):
        '''
        :param num: the wanted number of Schmidt weights.
        :param bond_head_site: the head site of the partition bond.
        :return: entanglement spectrum of self.
        '''
        if (bond_head_site == 0):
            return (-2) * np.log(self.v01[0: num])
        else:
            return (-2) * np.log(self.v10[0: num])

    def GetSchmidt(self, num=12, bond_head_site=0):
        '''
        :param num: the wanted number of Schmidt weights.
        :param bond_head_site: the head site of the partition bond.
        :return: Schmidt weights of self.
        '''
        if (bond_head_site == 0):
            return self.v01[0: num]
        else:
            return self.v10[0: num]

    def GetNorm(self, bond_head_site):
        '''
        :param bond_head_site: the head site of the wanted bond.
        :return: norm factor.
        '''
        state_tensor = self.ContractState(head_site=bond_head_site)
        state_norm = np.tensordot(
            state_tensor.conj(),
            state_tensor,
            axes=[[0, 1, 2, 3], [0, 1, 2, 3]],
        )
        # axes=[[1, 2], [0, 1]] means
        #   0---v---S----v----S---v---3
        #   :       |         |       :
        #   :       1         2       :
        #   :       :         :       :
        #   :       1         2       :
        #   :       |         |       :
        #   0---v---S----v----S---v---3
        return state_norm

    def MeasureBondE(self, bond_head_site):
        '''
        :param bond_head_site: the head site of the wanted bond.
        :return: bond energy.
        '''
        state_tensor = self.ContractState(head_site=bond_head_site)
        bond_E = np.tensordot(
            np.tensordot(
                state_tensor.conj(),
                self.H_pair,
                axes=[[1, 2], [0, 1]],
            ),
            state_tensor,
            axes=[[0, 2, 3, 1], [0, 1, 2, 3]],
        )
        # axes=[[1, 2], [0, 1]] means
        #   0---v---S----v----S---v---3
        #   :       |         |       :
        #   :       1         2       :
        #   :       :         :       :
        #   :       0         1       :
        #   :       |         |       :
        #   :       H---------H       :
        #   :       |         |       :
        #   :       2         3       :
        #   :       :         :       :
        #   :       1         2       :
        #   :       |         |       :
        #   0---v---S----v----S---v---3
        # axes=[[0, 2, 3, 1], [0, 1, 2, 3]] means
        #   0---v---S----v----S---v---1
        #   :       |         |       :
        #   :       |         |       :
        #   :       |         |       :
        #   :       |         |       :
        #   :       |         |       :
        #   :       H---------H       :
        #   :       |         |       :
        #   :       2         3       :
        #   :       :         :       :
        #   :       1         2       :
        #   :       |         |       :
        #   0---v---S----v----S---v---3
        # normalization, in fact it is already normalized if running correctly
        state_norm = self.GetNorm(bond_head_site=bond_head_site)
        bond_E /= state_norm
        return bond_E

    def MeasureAveE(self):
        '''
        :return: energy expectation averaged over bonds
        '''
        return 0.5 * (self.MeasureBondE(bond_head_site=0) + self.MeasureBondE(bond_head_site=1))


def GetPairHam_Heisenberg(J, dim_spin):
    '''
    :param J: coupling constant.
    :param dim_spin: Hilbert space dimension of single spin.
    :return: single two-site interaction term as a 4-order tensor J * S_i · S_j.
    '''
    Sz, Sp, Sm, Sx, Id = GenSpinOpr(which_opr='real', dim_spin=dim_spin)
    # H = J * \sum{ S_i * S_j }
    H_pair = np.kron(Sz, Sz) + 0.5 * J * (np.kron(Sp, Sm) + np.kron(Sm, Sp))
    #     0
    #    / \
    #   /   \
    #  H     H
    #   \   /
    #    \ /
    #     1
    H_pair = np.reshape(H_pair, [dim_spin, dim_spin, dim_spin, dim_spin])
    #  0     1
    #  |     |
    #  H-----H
    #  |     |
    #  2     3
    return H_pair


def GetPairHam_HaldanePhase(J=1, Uzz=0, Bx=0, Bz=0, R=0, dim_spin=3):
    '''
    :param J: coupling constant.
    :param Uzz: spin suppression constant along z direction.
    :param Bx: Zeeman field along x.
    :param Bz: Zeeman field along z.
    :param R: inversion symmetry breaking term constant.
    :param dim_spin: Hilbert space dimension of single spin.
    :return: single two-site interaction term as a 4-order tensor.
    '''
    Sz, Sp, Sm, Sx, Id = GenSpinOpr(which_opr='real', dim_spin=dim_spin)
    Sz2 = np.dot(Sz, Sz)
    R_term = np.kron(np.dot(Sz, Sx), Sx) - np.kron(Sx, np.dot(Sz, Sx))
    # H = J * \sum{ S_i * S_j } + Uzz * \sum{ (S^z_i)^2 } + Bx * \sum{ S^x_i }
    H_pair = np.kron(Sz, Sz) + 0.5 * J * (np.kron(Sp, Sm) + np.kron(Sm, Sp)) \
             + Uzz * np.kron(Sz2, Id) + Bx * np.kron(Sx, Id) + Bz * np.kron(Sz, Id) + R * R_term
    #     0
    #    / \
    #   /   \
    #  H     H
    #   \   /
    #    \ /
    #     1
    H_pair = np.reshape(H_pair, [dim_spin, dim_spin, dim_spin, dim_spin])
    #  0     1
    #  |     |
    #  H-----H
    #  |     |
    #  2     3
    return H_pair


def GetExpHamPair(H_pair, dt):
    '''
    :param H_pair: two-site operator as a 4-order tensor.
    :param dt: imaginary time evolution step.
    :return: a 4-order tensor of exp(- H_pair * dt).
    '''
    # exp(- H_pair * dt) = U * exp(- E_pair * dt) * UT
    d = H_pair.shape[0]
    # eigen-value decomposition
    E_pair, U = np.linalg.eigh(
        np.reshape(H_pair, [d * d, d * d])
    )
    # take the exponent of eigen values
    exp_Et = np.diag(np.exp(- E_pair * dt))
    # put them back together
    exp_Ht = np.reshape(
        np.dot(
            np.dot(
                U,
                exp_Et,
            ),
            U.conj().T,
        ),
        [d, d, d, d],
    )

    return exp_Ht