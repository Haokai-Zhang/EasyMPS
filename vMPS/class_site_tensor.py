#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <class_site_tensor.py> defines the class of tensors of single site.
'''

import numpy as np
import scipy.sparse.linalg

class site_tensor(object):
    '''
    class of state tensor combined with Hamiltonian mpo corresponding to a single site.
    '''

    # state tensor for a single site: 3-order tensor (2 virtual + 1 physical). Marked as "S" for "State".
    #             shape[1]                  (d)                                  (d)             (d)
    #               |                        |         for boundary sites:        |               |
    #   shape[0] ---S--- shape[2]      (D)---S---(D)                        (1)---S---(D)   (D)---S---(1)
    # we choose '012' instead of '021' for convenience consideration of reshape operation later.
    # the same for its Hermitian conjugate
    #   shape[0] ---S--- shape[2]      (D)---S---(D)                        (1)---S---(D)   (D)---S---(1)
    #               |                        |         for boundary sites:        |               |
    #             shape[1]                  (d)                                  (d)             (d)

    # the ASCII sketch of a series of site_tensor
    #      |   |   |   |   |   |   |   |   |
    #  e---S---S---S---S---S---S---S---S---S---e
    # 'e' is called an endpoint for programming convenience

    # operator tensor for a single site: 4-order tensor (2 virtual + 2 physical). Marked as "H" for "Hamiltonian".
    #             shape[2]                  (d)                                   (d)             (d)
    #               |                        |          for boundary sites:        |               |
    #   shape[0] ---H--- shape[1]      (D)---H---(D)                         (1)---H---(D)   (D)---H---(1)
    #               |                        |                                     |               |
    #             shape[3]                  (d)                                   (d)             (d)
    # (D),(d) represents bond dimension (the number of summation terms).

    # initialize
    def __init__(self, D_bond_l, D_bond_r, mpo):
        '''
        :param D_bond_l: left bond dimension of site_tensor.
        :param D_bond_r: right bond dimension of site_tensor.
        :param mpo: Hamiltonian mpo corresponding to site_tensor.
        :return: initialized site_tensor.
        '''

        # site number
        self.num = 0

        # pointer to the left/right neighbor initialized as None
        self.neighbor_l = None
        self.neighbor_r = None

        if (D_bond_l > 0) and (D_bond_r > 0):
            # if it is not an endpoint
            self.is_end = False
            # random initialization
            self.tensor_data = np.random.rand(D_bond_l, mpo.shape[2], D_bond_r)
            self.mpo = mpo
            # save cache to accelerate calculations at the cost of memory
            self.M_cache = None
            self.L_cache = None
            self.R_cache = None
            self.sv = None
            # E = <psi|H|psi> =
            #   (1)---S---S---S---S---S---S---S---S---(1)
            #         |   |   |   |   |   |   |   |
            #   (1)---H---H---H---H---H---H---H---H---(1)
            #         |   |   |   |   |   |   |   |
            #   (1)---S---S---S---S---S---S---S---S---(1)
            # which has 6 legs of dimension = 1, i.e. a scalar in fact.
            # denote tensor block L, M, R as
            #  ---S---S---S---  ---S---  ---S---S---S---S---                   ---|'''|---  ---S---  ---|'''|---
            #     |   |   |        |        |   |   |   |                         |   |        |        |   |
            #  ---H---H---H---  ---H---  ---H---H---H---H---   ---briefly--->  ---|   |---  ---H---  ---|   |---
            #     |   |   |        |        |   |   |   |                         |   |        |        |   |
            #  ---S---S---S---  ---S---  ---S---S---S---S---                   ---|...|---  ---S---  ---|...|---
        else:
            # if it is an endpoint
            self.is_end = True

    # give the dimension of each bond a special name
    # to facilitate calling
    @property
    def D_bond_l(self):
        return self.data_protected.shape[0]
    @property
    def d_phys(self):
        return self.data_protected.shape[1]
    @property
    def D_bond_r(self):
        return self.data_protected.shape[2]

    # define a function of inserting a new site nearby
    def InsertNewSite(self, new_site, side_lr):
        '''
        :param new_site: initialized new site_tensor
        :param side_lr: {"l", "r"} determine the left/right side to insert by.
        :return: None. Insert 'new_site' on the side of 'self' by changing the neighboring relation.
        '''
        if (side_lr == "l"):
            # record the original left neighbor
            l_temp = self.neighbor_l
            #  l---self---r  =>  l-->self---r
            #                        /
            #                    new_site
            self.neighbor_l = new_site
            new_site.neighbor_r = self
            #  l-->self---r  =>  l   self---r
            #      /              \  /
            #  new_site          new_site
            l_temp.neighbor_r = new_site
            new_site.neighbor_l = l_temp
        elif (side_lr == "r"):
            # record the original tight neighbor
            r_temp = self.neighbor_r
            #  l---self---r  =>  l---self<--r
            #                           \
            #                         new_site
            self.neighbor_r = new_site
            new_site.neighbor_l = self
            #  l---self<--r  =>  l---self   r
            #         \                 \  /
            #       new_site          new_site
            r_temp.neighbor_l = new_site
            new_site.neighbor_r = r_temp
        else:
            print("----Error: (side_lr = l or r) violated !----")
            exit(1)



    # protection against writing in
    # clear cache as long as self.tensor_data changes
    @property
    def data_protected(self):
        return self.tensor_data
    @data_protected.setter
    def data_protected(self, new_state_tensor):
        self.tensor_data= new_state_tensor
        # clear cache involving self.data_protected
        self.ClearCacheLR()

    # define functions of cache clear
    def ClearCacheLR(self):
        '''
        :return: None. Clear cache involving self.data_protected if it changes.
        '''
        self.M_cache = None
        self.neighbor_r.ClearCacheL()
        self.neighbor_l.ClearCacheR()

    def ClearCacheL(self):
        '''
        :return: None. Clear all L_cache for sites on the right by recursion.
        '''
        if self.is_end:
            # stop if meet the endpoint
            return
        self.L_cache = None
        self.neighbor_r.ClearCacheL()

    def ClearCacheR(self):
        '''
        :return: None. Clear all R_cache for sites on the left by recursion.
        '''
        if self.is_end:
            # stop if meet the endpoint
            return
        self.R_cache = None
        self.neighbor_l.ClearCacheR()




    # define singular value decomposition functions for left/right bond
    def SvdBondL(self):
        '''
        :return: None. Svd for the left bond of self to make self satisfy right canonical condition.
        '''

        # the sketch of the left bond of 'self'.
        #       |   |   |   |   |            |     |   |   |
        #    ---S---S---S---S---l--<svd>--(self)---r---S---S---

        # return if meet the endpoint
        if self.is_end:
            return

        # implement matrix form svd after reshaping
        U, S, VT = np.linalg.svd(self.data_protected.reshape(self.D_bond_l,
                                                             self.d_phys * self.D_bond_r),
                                 full_matrices=False,
                                 )
        # full_matrices=False means the reduced svd,
        # i.e. U=M*K VT=K*N K=min(M,N) instead of M*M and N*N

        # recognize VT as the new data of 'self'
        # so that right canonical condition can be satisfied
        #    .---S---     .---
        #    |   |     =  |     =  1_{D*D}
        #    '---S---     '---
        self.data_protected = VT.reshape((-1, self.d_phys, self.D_bond_r))
        # parameter '-1' means that it will be automatically calculated
        # according to the other parameters

        # recognize l_old*U*S as the new data of the left neighbor of 'self'
        if (not self.neighbor_l.is_end):
            self.neighbor_l.data_protected = np.tensordot(self.neighbor_l.data_protected,
                                                          np.dot(U, np.diag(S)),
                                                          axes=[2, 0],
                                                          )
        # 'axes=[2, 0]' means
        #        1
        #        |
        #   0 ---l--- 2 <----> 0 ---U*S--- 1
        # 'numpy.diag()' turns a vector into a diagonal matrix

    def SvdBondR(self):
        '''
        :return: None. Svd for the right bond of self to make self satisfy left canonical condition.
        '''

        # the sketch of the right bond of 'self'.
        #       |   |   |   |      |           |   |   |   |
        #    ---S---S---S---l---(self)--<svd>--r---S---S---S---

        # return if meet the endpoint
        if self.is_end:
            return

        # implement matrix form svd after reshaping
        U, S, VT = np.linalg.svd(self.data_protected.reshape(self.D_bond_l * self.d_phys,
                                                             self.D_bond_r),
                                 full_matrices=False,
                                 )
        # full_matrices=False means the reduced svd,
        # i.e. U=M*K VT=K*N K=min(M,N) instead of M*M and N*N

        # recognize U as the new data of 'self'
        # so that left canonical condition can be satisfied
        #    ---S---.     ---.
        #       |   |  =     |  =  1_{D*D}
        #    ---S---'     ---'
        self.data_protected = U.reshape((self.D_bond_l, self.d_phys, -1))
        # parameter '-1' means that it will be automatically calculated
        # according to the other parameters

        # recognize S*VT*r_old as the new data of the right neighbor of 'self'
        if (not self.neighbor_r.is_end):
            self.neighbor_r.data_protected = np.tensordot(np.dot(np.diag(S), VT),
                                                          self.neighbor_r.data_protected,
                                                          axes=[1, 0],
                                                          )
        # 'axes=[1, 0]' means
        #                              1
        #                              |
        #   0 ---S*VT--- 1 <----> 0 ---r--- 2
        # note that '1' in '--- 1' means shape[1]
        # 'numpy.diag()' turns a vector into a diagonal matrix



    # define functions of calculating block L, M, R, which can be sketched as
    #  ---S---S---S---  ---S---  ---S---S---S---S---                   ---|'''|---  ---S---  ---|'''|---
    #     |   |   |        |        |   |   |   |                         |   |        |        |   |
    #  ---H---H---H---  ---H---  ---H---H---H---H---   ---briefly--->  ---|   |---  ---H---  ---|   |---
    #     |   |   |        |        |   |   |   |                         |   |        |        |   |
    #  ---S---S---S---  ---S---  ---S---S---S---S---                   ---|...|---  ---S---  ---|...|---
    # note that the general contraction index rules: "number them in original order"
    # and "the more before, the larger".
    def GetBlockM(self):
        '''
        :return: block M of self, i.e.an sandwich-shaped 'S-H-S' block.
        '''
        if (self.M_cache is not None):
            # return cache if calculated before and self.data_protected has never changed
            return self.M_cache
        else:
            # re-calculate if no cache
            S_H = np.tensordot(self.data_protected.conj(),
                               self.mpo,
                               axes=[1, 2])
            S_H_S = np.tensordot(S_H,
                                 self.data_protected,
                                 axes=[4, 1])
            # 'axes=[1, 2]' and 'axes=[4, 1]' mean
            #   0 ---S--- 2                  0 ---S--- 1                  0 ---S--- 1
            #        |                            |                            |
            #        1        contraction         |                            |
            #                 ===========>        |                            |
            #        2                            |                            |
            #        |                            |                            |
            #   0 ---H--- 1                  2 ---H--- 3                  2 ---H--- 3
            #        |                            |                            |
            #        3                            4        contraction         |
            #                                              ===========>        |
            #        1                            1                            |
            #        |                            |                            |
            #   0 ---S--- 2                  0 ---S--- 2                  4 ---S--- 5
            self.M_cache = S_H_S
            return self.M_cache

    def GetBlockL(self):
        '''
        :return: block L of self (recursion). self.GetBlockL() excludes 'self'.
        '''
        if self.neighbor_l.is_end:
            # return scalar 1 with 6-pseudo-index if meet endpoint
            #     (1) (1)
            #       \ /
            #  (1)---1---(1)
            #       / \
            #     (1) (1)
            # note that A---(1)---B actually means the tensor product of A and B with no contraction
            return np.ones([1, 1, 1, 1, 1, 1])
        if (self.L_cache is not None):
            # return cache if calculated before and self.data_protected has never changed
            return self.L_cache
        else:
            # re-calculate if no cache
            self.L_cache = np.tensordot(self.neighbor_l.GetBlockL(),
                                        self.neighbor_l.GetBlockM(),
                                        axes=[[1, 3, 5], [0, 2, 4]],
                                        ).transpose([0, 3, 1, 4, 2, 5])
            # axes=[[1, 3, 5], [0, 2, 4]] means
            #   0 ---S--- 1  <--->  0 ---S--- 1                  0 ---S---S--- 3
            #        |                   |        contraction         |   |
            #   2 ---H--- 3  <--->  2 ---H--- 3   ===========>   1 ---H---H--- 4
            #        |                   |                            |   |
            #   4 ---S--- 5  <--->  4 ---S--- 5                  2 ---S---S--- 5
            # transpose([0, 3, 1, 4, 2, 5]) means
            #   0 ---|'''|--- 3                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|   |--- 4   =========>   2 ---|   |--- 3
            #        |   |                          |   |
            #   2 ---|   |--- 5                4 ---|...|--- 5
            return self.L_cache

    def GetBlockR(self):
        '''
        :return: block R of self (recursion). self.GetBlockR() excludes 'self'.
        '''
        if self.neighbor_r.is_end:
            # return scalar 1 with 6-pseudo-index if meet endpoint, i.e.
            #     (1) (1)
            #       \ /
            #  (1)---1---(1)
            #       / \
            #     (1) (1)
            # note that A---(1)---B actually means the tensor product of A and B with no contraction
            return np.ones([1, 1, 1, 1, 1, 1])
        if (self.R_cache is not None):
            # return cache if calculated before and self.data_protected has never changed
            return self.R_cache
        else:
            self.R_cache = np.tensordot(self.neighbor_r.GetBlockM(),
                                        self.neighbor_r.GetBlockR(),
                                        axes=[[1, 3, 5], [0, 2, 4]],
                                        ).transpose([0, 3, 1, 4, 2, 5])
            # axes=[[1, 3, 5], [0, 2, 4]] means
            #   0 ---S--- 1  <--->  0 ---S--- 1                  0 ---S---S--- 3
            #        |                   |        contraction         |   |
            #   2 ---H--- 3  <--->  2 ---H--- 3   ===========>   1 ---H---H--- 4
            #        |                   |                            |   |
            #   4 ---S--- 5  <--->  4 ---S--- 5                  2 ---S---S--- 5
            # transpose([0, 3, 1, 4, 2, 5]) means
            #   0 ---|'''|--- 3                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|   |--- 4   =========>   2 ---|   |--- 3
            #        |   |                          |   |
            #   2 ---|   |--- 5                4 ---|...|--- 5
            return self.R_cache


    # define a function of calculating the effective Hamiltonian for the 2 sites to be updated
    def GetHeff2Sites(self, update_with):
        '''
        :param update_with: {"l", "r"} determine the left/right side to update together with.
        :return: the effective Hamiltonian "L-mpo-mpo-R" corresponding to the 2 sites to be updated.
        '''

        # the sketch of the effective Hamiltonian H_eff
        #   |'''|---              ---|'''|
        #   |   |      |      |      |   |
        #   |   |------H------H------|   |
        #   |   |      |      |      |   |
        #   |...|---              ---|...|
        # as a 8-order tensor

        # choose the site_tensor of which side to update together with
        if (update_with == "l"):
            s0 = self.neighbor_l
            s1 = self
        elif (update_with == "r"):
            s0 = self
            s1 = self.neighbor_r
        else:
            print("----Error: (update_with = l or r) violated !----")
            exit(1)
        # the sketch of s0 and s1
        #       |   |   |   |   |    |    |   |
        #    ---S---S---S---S---s0---s1---S---S---

        # contract (L---mpo---mpo---R) to get H_eff
        H_eff_order_14 = np.tensordot(
                                    np.tensordot(
                                        np.tensordot(s0.GetBlockL(),
                                                     s0.mpo,
                                                     axes=[3, 0],
                                                     ),
                                        s1.mpo,
                                        axes=[5, 0],
                                    ),
                                    s1.GetBlockR(),
                                    axes=[7, 2],
                                )
        # axes=[3, 0] means (L---mpo)
        #  0 ---|'''|--- 1              2                     0 ---|'''|--- 1     6
        #       |   |                   |       contraction        |   |          |
        #  2 ---|   |--- 3  <--->  0 ---H--- 1  ===========>  2 ---|   |----------H--- 5
        #       |   |                   |                          |   |          |
        #  4 ---|...|--- 5              3                     3 ---|...|--- 4     7
        # axes=[5, 0] means (L-mpo---mpo)
        #  0 ---|'''|--- 1    6                   2                     0 ---|'''|--- 1     5      8
        #       |   |         |                   |       contraction        |   |          |      |
        #  2 ---|   |---------H--- 5  <--->  0 ---H--- 1  ===========>  2 ---|   |----------H------H--- 7
        #       |   |         |                   |                          |   |          |      |
        #  3 ---|...|--- 4    7                   3                     3 ---|...|--- 4     6      9
        # axes=[7, 2] means (L-mpo-mpo---R)
        #  0 ---|'''|--- 1     5      8              0 ---|'''|--- 1
        #       |   |          |      |                   |   |
        #  2 ---|   |----------H------H--- 7  <--->  2 ---|   |--- 3
        #       |   |          |      |                   |   |
        #  3 ---|...|--- 4     6      9              4 ---|...|--- 5
        #
        #                  0 ---|'''|--- 1     5      7     9 ---|'''|--- 10
        #   contraction         |   |          |      |          |   |
        #   ===========>   2 ---|   |----------H------H----------|   |--- 11
        #                       |   |          |      |          |   |
        #                  3 ---|...|--- 4     6      8    12 ---|...|--- 13

        # reduce the 6-pseudo-index and transpose for the next step
        H_eff_order_8 = H_eff_order_14.reshape(
                                    s0.D_bond_l,
                                    s0.D_bond_l,
                                    s0.d_phys,
                                    s0.d_phys,
                                    s1.d_phys,
                                    s1.d_phys,
                                    s1.D_bond_r,
                                    s1.D_bond_r,
                                ).transpose([0, 2, 4, 6, 1, 3, 5, 7])
        # reshape(D, D, d, d, d, d, D, D) means
        #                |'''|--- 0     2      4     6 ---|'''|
        #    reduce      |   |          |      |          |   |
        #    ======>     |   |----------H------H----------|   |
        #                |   |          |      |          |   |
        #                |...|--- 1     3      5     7 ---|...|
        # transpose([0, 2, 4, 6, 1, 3, 5, 7]) means
        #                |'''|--- 0     1      2     3 ---|'''|
        #   transpose    |   |          |      |          |   |
        #   =========>   |   |----------H------H----------|   |
        #                |   |          |      |          |   |
        #                |...|--- 4     5      6     7 ---|...|

        # reshape H_eff to a square matrix
        dim_H_eff_matrix = s0.D_bond_l * s0.d_phys * s1.d_phys * s1.D_bond_r
        H_eff_order_2 = H_eff_order_8.reshape(dim_H_eff_matrix, dim_H_eff_matrix)

        return H_eff_order_2

    # define a function of 2-site update
    def Update2Sites(self, towards, update_with, if_print=False):
        '''
        :param towards: {"l", "r"} determine the sweep direction.
        :param update_with: {"l", "r"} determine the left/right side to update together with.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :return: the lowest energy after this update.
        '''

        # diagonalize H_eff to get the ground state
        # the eigen-vector is normalized by default of 'eigsh()'
        E_complex, V_complex = scipy.sparse.linalg.eigsh(self.GetHeff2Sites(update_with),
                                                         1,
                                                         which="SA",
                                                         )

        # take the real part to avoid the infinitesimal imaginary part
        E_real = E_complex[0].real

        if (update_with == "l"):
            s0 = self.neighbor_l
            s1 = self
        elif (update_with == "r"):
            s0 = self
            s1 = self.neighbor_r
        else:
            print("----Error: (update_with = l or r) violated !----")
            exit(1)
        # the sketch of s0 and s1
        #       |   |   |   |   |    |    |   |
        #    ---S---S---S---S---s0---s1---S---S---

        # reshape the eigen-vector to a square matrix for svd
        state_tensor_2_sites = V_complex.reshape(s0.D_bond_l * s0.d_phys,
                                                 s1.d_phys * s1.D_bond_r,
                                                 )

        # svd the 4-order tensor to two 3-order tensors
        U, S, VT = np.linalg.svd(state_tensor_2_sites, full_matrices=False)

        # truncate
        # the number of singular values, i.e.the bond dimension between U&S*VT / U*S&VT
        # might be larger than D due to D*d>D
        D = s0.D_bond_r
        # evaluate the truncation error
        trunc_err = 0
        if (len(S) > D):
            trunc_err = S[D:].sum()
            U = U[:, 0:D]
            S = S[0:D]
            VT = VT[0:D, :]
        # record the singular values
        s0.sv = S

        if (towards == "l"):
            s0.data_protected = np.dot(U,
                                       np.diag(S),
                                       ).reshape((s0.D_bond_l,
                                                  s0.d_phys,
                                                  -1,
                                                  ))
            s1.data_protected = VT.reshape((-1,
                                            s1.d_phys,
                                            s1.D_bond_r,
                                            ))
            # recover the left/right canonical condition
            s0.SvdBondL()
        elif (towards == "r"):
            s0.data_protected = U.reshape((s0.D_bond_l,
                                           s0.d_phys,
                                           -1,
                                           ))
            s1.data_protected = np.dot(np.diag(S),
                                       VT,
                                       ).reshape((-1,
                                                  s1.d_phys,
                                                  s1.D_bond_r,
                                                  ))
            s1.SvdBondR()
        else:
            print("----Error: (towards = l or r) violated !----")
            exit(1)

        # print the procedure of sweep
        if (if_print == True):
            print('site', end='')
            site_str = '(%d,%d)' % (s0.num,s1.num)
            print('%7s' % site_str, end=' ')
            print('E=%.14f' % E_real, end=' ')
            print('D=%3d' % s0.D_bond_r, end=' ')
            print('TruncErr=', format(trunc_err, '.1e'))

        return E_real
    
    def GetHeff1Site(self):
        '''
        :return: the effective Hamiltonian "L-mpo-R" corresponding to the site to be updated.
        '''

        # the sketch of the effective Hamiltonian H_eff
        #   |'''|---       ---|'''|
        #   |   |      |      |   |
        #   |   |------H------|   |
        #   |   |      |      |   |
        #   |...|---       ---|...|
        # as a 6-order tensor

        # contract (L---mpo---R) to get H_eff
        H_eff_order_12 = np.tensordot(
                                    np.tensordot(
                                        self.GetBlockL(),
                                        self.mpo,
                                        axes=[3, 0],
                                    ),
                                    self.GetBlockR(),
                                axes=[5, 2],
                                )
        # axes=[3, 0] means (L---mpo)
        #  0 ---|'''|--- 1              2                     0 ---|'''|--- 1     6
        #       |   |                   |       contraction        |   |          |
        #  2 ---|   |--- 3  <--->  0 ---H--- 1  ===========>  2 ---|   |----------H--- 5
        #       |   |                   |                          |   |          |
        #  4 ---|...|--- 5              3                     3 ---|...|--- 4     7
        # axes=[5, 2] means (L-mpo-mpo---R)
        #  0 ---|'''|--- 1     6              0 ---|'''|--- 1
        #       |   |          |                   |   |
        #  2 ---|   |----------H--- 5  <--->  2 ---|   |--- 3
        #       |   |          |                   |   |
        #  3 ---|...|--- 4     7              4 ---|...|--- 5
        #
        #                  0 ---|'''|--- 1     5     7 ---|'''|--- 8
        #   contraction         |   |          |          |   |
        #   ===========>   2 ---|   |----------H----------|   |--- 9
        #                       |   |          |          |   |
        #                  3 ---|...|--- 4     6    10 ---|...|--- 11

        # reduce the 6-pseudo-index and transpose for the next step
        H_eff_order_6 = H_eff_order_12.reshape(self.D_bond_l,
                                               self.D_bond_l,
                                               self.d_phys,
                                               self.d_phys,
                                               self.D_bond_r,
                                               self.D_bond_r,
                                            ).transpose([0, 2, 4, 1, 3, 5])
        # reshape(D, D, d, d, D, D) means
        #                |'''|--- 0     2     4 ---|'''|
        #    reduce      |   |          |          |   |
        #    ======>     |   |----------H----------|   |
        #                |   |          |          |   |
        #                |...|--- 1     3     5 ---|...|
        # transpose([0, 2, 4, 1, 3, 5]) means
        #                |'''|--- 0     1     2 ---|'''|
        #   transpose    |   |          |          |   |
        #   =========>   |   |----------H----------|   |
        #                |   |          |          |   |
        #                |...|--- 3     4     5 ---|...|

        # reshape H_eff to a square matrix
        dim_Heff_matrix = self.D_bond_l * self.d_phys * self.D_bond_r
        H_eff_order_2 = H_eff_order_6.reshape(dim_Heff_matrix, dim_Heff_matrix)

        return H_eff_order_2

    # define a function of 2-site update
    def Update1Site(self, towards, if_print=False):
        '''
        :param towards: {"l", "r"} determine the sweep direction.
        :param if_print: {True, False} determine whether to print the sweep procedure.
        :return: the lowest energy after this update.
        '''

        # diagonalize H_eff to get the ground state
        # the eigen-vector is normalized by default of 'eigsh()'
        E_complex, V_complex = scipy.sparse.linalg.eigsh(self.GetHeff1Site(),
                                                         1,
                                                         which="SA",
                                                         )

        # take the real part to avoid the infinitesimal imaginary part
        E_real = E_complex[0].real

        # reshape the vector back to a 3-order state tensor
        self.data_protected = V_complex.reshape(self.D_bond_l,
                                                self.d_phys,
                                                self.D_bond_r,
                                                )

        # recover the left/right canonical condition
        if (towards == "l"):
            self.SvdBondL()
        elif (towards == "r"):
            self.SvdBondR()
        else:
            print("----Error: (towards = l or r) violated !----")
            exit(1)

        # print the procedure of sweep
        if (if_print == True):
            print('site', end='')
            site_str = '%d' % self.num
            print('%3s' % site_str, end=' ')
            print('E=%.14f' % E_real, end=' ')
            if (towards == "l"):
                print('D=%d' % self.D_bond_l)
            if (towards == "r"):
                print('D=%d' % self.D_bond_r)

        return E_real

    def GetNormM(self):
        '''
        :return: the norm tensor of self.
        '''
        if (not self.is_end):
            S_S = np.tensordot(self.data_protected.conj(),
                               self.data_protected,
                               axes=[1, 1])
            # 'axes=[1, 1]' means
            #   0 ---S--- 2                  0 ---S--- 1
            #        |                            |
            #        1                            |
            #               --contraction-->      |
            #        1                            |
            #        |                            |
            #   0 ---S--- 2                  2 ---S--- 3
            return S_S
        else:
            print("----Error: endpoint has no norm !----")
            exit(1)

    def GetNormL(self):
        '''
        :return: the norm tensor of the left part of self by recursion. self.GetNormL() includes self.
        '''
        if self.is_end:
            # return scalar 1 with 4-pseudo-index if meet endpoint
            return np.ones([1, 1, 1, 1])
        else:
            norm_L = np.tensordot(self.neighbor_l.GetNormL(),
                                  self.GetNormM(),
                                  axes=[[1, 3], [0, 2]],
                                  ).transpose([0, 2, 1, 3])
            # axes=[[1, 3], [0, 2]] means
            #   0 ---S--- 1  <--->  0 ---S--- 1                  0 ---S---S--- 2
            #        |                   |        contraction         |   |
            #   2 ---S--- 3  <--->  2 ---S--- 3   ===========>   1 ---S---S--- 3
            # transpose([0, 2, 1, 3]) means
            #   0 ---|'''|--- 2                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|...|--- 3   =========>   2 ---|...|--- 3
            return norm_L

    def GetNormR(self):
        '''
        :return: the norm tensor of the right part of self by recursion. self.GetNormR() includes self.
        '''
        if self.is_end:
            # return scalar 1 with 4-pseudo-index if meet endpoint
            return np.ones([1, 1, 1, 1])
        else:
            norm_R = np.tensordot(self.GetNormM(),
                                  self.neighbor_r.GetNormR(),
                                  axes=[[1, 3], [0, 2]],
                                  ).transpose([0, 2, 1, 3])
            # axes=[[1, 3], [0, 2]] means
            #   0 ---S--- 1  <--->  0 ---S--- 1                  0 ---S---S--- 2
            #        |                   |        contraction         |   |
            #   2 ---S--- 3  <--->  2 ---S--- 3   ===========>   1 ---S---S--- 3
            # transpose([0, 2, 1, 3]) means
            #   0 ---|'''|--- 2                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|...|--- 3   =========>   2 ---|...|--- 3
            return norm_R

    def GetMeasureM(self, opr):
        '''
        :param opr: local operator to be measured
        :return: a sandwich-shaped tensor obtained by contracting S-opr-S
        '''
        if (len(opr.shape) != 2):
            print("----Error from GetMeasureM(): len(opr.shape) == 2 violated !----")
            exit(1)
        if (not self.is_end):
            S_O = np.tensordot(self.data_protected.conj(),
                               opr,
                               axes=[1, 0])
            S_O_S = np.tensordot(S_O,
                                 self.data_protected,
                                 axes=[2, 1])
            # 'axes=[1, 0]' and 'axes=[2, 1]' mean
            #   0 ---S--- 2                  0 ---S--- 1                  0 ---S--- 1
            #        |                            |                            |
            #        1        contraction         |                            |
            #                 ===========>        |                            |
            #        0                            |                            |
            #        |                            |                            |
            #        O                            O                            O
            #        |                            |                            |
            #        1                            2        contraction         |
            #                                              ===========>        |
            #        1                            1                            |
            #        |                            |                            |
            #   0 ---S--- 2                  0 ---S--- 2                  2 ---S--- 3
            # where '0---O---1' denotes the local operator 'opr'
            return S_O_S
        else:
            print("----Error from GetMeasureM(): endpoint !----")
            exit(1)

    def GetMeasureL(self, list_opr, list_site):
        '''
        :param list_opr: list of local operators to be measured on the left of self (including self)
        :param list_site: list of sites corresponding to these operators on the left of self (including self)
        :return: a tensor obtained by contracting the operators with state tensors on the left of self (including self)
        '''
        if self.is_end:
            # return scalar 1 with 4-pseudo-index if meet endpoint
            return np.ones([1, 1, 1, 1])
        else:
            if (len(list_site) == 0):
                return self.GetNormL()
            elif (self.num == list_site[-1]):
                # if there is an operator corresponding to self
                # take it out from list_site and list_opr and measure it
                list_site.pop(-1)
                opr = list_opr.pop(-1)
                # contract the measured result of self with its left (recursion)
                Measure_L = np.tensordot(self.neighbor_l.GetMeasureL(list_opr, list_site),
                                         self.GetMeasureM(opr=opr),
                                         axes=[[1, 3], [0, 2]],
                                         ).transpose([0, 2, 1, 3])
            else:
                # if there is no operator corresponding to self but yes for the left of self
                # take the norm of this site directly and measure the next
                Measure_L = np.tensordot(self.neighbor_l.GetMeasureL(list_opr, list_site),
                                         self.GetNormM(),
                                         axes=[[1, 3], [0, 2]],
                                         ).transpose([0, 2, 1, 3])
            # axes=[[1, 3], [0, 2]] means
            #   0 ---|'''|--- 1  <--->  0 ---|'''|--- 1                  0 ---|'''|--- 2
            #        |   |                   |   |        contraction         |   |
            #   2 ---|...|--- 3  <--->  2 ---|...|--- 3   ===========>   1 ---|...|--- 3
            # transpose([0, 2, 1, 3]) means
            #   0 ---|'''|--- 2                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|...|--- 3   =========>   2 ---|...|--- 3
            return Measure_L

    def GetMeasureR(self, list_opr, list_site):
        '''
        :param list_opr: list of local operators to be measured on the right of self (including self)
        :param list_site: list of sites corresponding to these operators on the right of self (including self)
        :return: a tensor obtained by contracting the operators with state tensors on the right of self (including self)
        '''
        if self.is_end:
            # return scalar 1 with 4-pseudo-index if meet endpoint
            return np.ones([1, 1, 1, 1])
        else:
            if (len(list_site) == 0):
                # if there is no more operators to measure on the right of self (including self)
                # return the norm directly
                return self.GetNormR()
            elif (self.num == list_site[0]):
                # if there is an operator corresponding to self
                # take it out from list_site and list_opr and measure it
                list_site.pop(0)
                opr = list_opr.pop(0)
                # contract the measured result of self with its right (recursion)
                Measure_R = np.tensordot(self.GetMeasureM(opr=opr),
                                         self.neighbor_r.GetMeasureR(list_opr, list_site),
                                         axes=[[1, 3], [0, 2]],
                                         ).transpose([0, 2, 1, 3])
            else:
                # if there is no operator corresponding to self but yes for the right of self
                # take the norm of this site directly and measure the next
                Measure_R = np.tensordot(self.GetNormM(),
                                         self.neighbor_r.GetMeasureR(list_opr, list_site),
                                         axes=[[1, 3], [0, 2]],
                                         ).transpose([0, 2, 1, 3])
            # axes=[[1, 3], [0, 2]] means
            #   0 ---|'''|--- 1  <--->  0 ---|'''|--- 1                  0 ---|'''|--- 2
            #        |   |                   |   |        contraction         |   |
            #   2 ---|...|--- 3  <--->  2 ---|...|--- 3   ===========>   1 ---|...|--- 3
            # transpose([0, 2, 1, 3]) means
            #   0 ---|'''|--- 2                0 ---|'''|--- 1
            #        |   |        transpose         |   |
            #   1 ---|...|--- 3   =========>   2 ---|...|--- 3
            return Measure_R

    # Supplemental comments:
    # The function of <np.array.reshape> is to rearrange the array in the same order as before
    # for example: (4, 2, 5) ---> (4*2, 5)
    # [[[1,2,3,4,5],[6,7,8,9,10]],                       [[1,2,3,4,5],[6,7,8,9,10],
    # [[11,12,13,14,15],[16,17,18,19,20]],     --->      [11,12,13,14,15],[16,17,18,19,20],
    # [[21,22,23,24,25],[26,27,28,29,30]],               [21,22,23,24,25],[26,27,28,29,30],
    # [[31,32,33,34,35],[36,37,38,39,40]]]               [31,32,33,34,35],[36,37,38,39,40]]
    # which is equivalent with taking the second bracket away
    # (4, 2, 5) ---> (4, 2*5)
    # [[[1,2,3,4,5],[6,7,8,9,10]],                       [[1,2,3,4,5,6,7,8,9,10],
    # [[11,12,13,14,15],[16,17,18,19,20]],     --->      [11,12,13,14,15,16,17,18,19,20],
    # [[21,22,23,24,25],[26,27,28,29,30]],               [21,22,23,24,25,26,27,28,29,30],
    # [[31,32,33,34,35],[36,37,38,39,40]]]               [31,32,33,34,35,36,37,38,39,40]]
    # which is equivalent with taking the third bracket away