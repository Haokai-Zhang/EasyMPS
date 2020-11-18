#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. Run <test_excited.py> to calculate the first excited state.
'''

from Test.model_inout import GetParams, PrintE
from EasyOpr.get_mpo import *
from vMPS.class_mps import mps
from EasyPlot.plot_xy import *
import numpy as np

def GetExcitedHmpo(list_mpo_H, mps0, E0):
    '''
    :param list_mpo_H: mpo list of the original Hamiltonian.
    :param mps0: mps class of the ground state.
    :param E0: ground state energy.
    :return: the variation mpo for the first excited state.
    '''
    # H' = H - E0 * | 0 > < 0 |

    # extract mps data
    list_mps_data = []
    the_site_tensor = mps0.end_l.neighbor_r
    while (not the_site_tensor.is_end):
        list_mps_data += [the_site_tensor.data_protected]
        the_site_tensor = the_site_tensor.neighbor_r

    # construct the ground state projector mpo by outer product for each site
    # note that this projector has bond dimension of D^2
    # which brings us contraction process of D^6 (ground state only D^4)
    list_mpo_projector = []
    for mps_data in list_mps_data:
        # append a pseudo-index for outer product
        new_shape = mps_data.shape + (1,)
        mps_data = np.reshape(mps_data, new_shape)
        #        1                           1
        #        |          pseudo           |
        #   0 ---S--- 2   ==========>   0 ---S--- 2
        #                    index           |
        #                                    3

        # outer product
        outer_prod = np.tensordot(mps_data,
                                  mps_data.conj(),
                                  axes=[3, 3],
                                  )
        # axes=[3, 3] means
        #        1                            1
        #        |                            |
        #   0 ---S--- 2                  0 ---S--- 2
        #        |                            |
        #        3        contraction         |
        #                 ===========>        |
        #        3                            |
        #        |                            |
        #   0 ---S--- 2                  3 ---S--- 5
        #        |                            |
        #        1                            4

        # combine the virtual dimensions of the projector
        # to make it consistent with the Hamiltonian mpo
        mpo_projector = CombineIndex(outer_prod,
                                     list_index_to_combine=[[0, 3],
                                                            [2, 5],
                                                            [1,],
                                                            [4,],
                                                            ])
        # list_combine_pair=[[0,3],[2,5],[1,],[4,]] means
        #        1                             2
        #        |                             |
        #   0 ---S--- 2                     ---S---
        #        |                         /   |   \                   2
        #        |        combine         /    |    \                  |
        #        |        =======>    0 --     |     -- 1   ===>  0 ---P--- 1
        #        |                        \    |    /                  |
        #        |                         \   |   /                   3
        #   3 ---S--- 5                     ---S---
        #        |                             |
        #        4                             3
        list_mpo_projector.append(mpo_projector)


    # projector multiplied by scalar (-E0)
    list_mpo_projector[0] = np.tensordot(np.asarray([[-E0]]),
                                         list_mpo_projector[0],
                                         axes=[1, 0],
                                         )
    # axes=[1, 0] means
    #                           2
    #                           |
    # 0 ---E0--- 1  <--->  0 ---P--- 1
    #                           |
    #                           3

    # construct the variation mpo for the first excited state by
    # adding the mpo of Hamiltonian and projector on the virtual bond,
    # H' = H - E0 * | 0 > < 0 |
    # where the additivity is guaranteed by the same physical bond dimension
    list_mpo_H_excited = []
    for site_index in range(0, len(list_mpo_H)):
        list_mpo_H_excited.append(ExtendAdd(list_mpo_H[site_index], list_mpo_projector[site_index]))

    return list_mpo_H_excited

def ExtendAdd(array_A, array_B):
    '''
    :param array_A: array to add
    :param array_B: array to add
    :return: an array obtained by add array_A and array_B on the indices with the same dimension.
    '''
    # example: ExtendAdd([[1,2],[3,4]], [[5,6],[7,8],[9,0]]) = [[1,2],[3,4],[5,6],[7,8],[9,0]]
    # example: ExtendAdd([[1,2],[3,4]], [[5,6,7],[8,9,0]]) = [[1,2,5,6,7],[3,4,8,9,0]]
    # which is exactly the same as np.concatenate()
    # example: ExtendAdd([[1,],[2,]], [[3,4]]) = [[1,0,0],[2,0,0],[0,3,4]]
    # example: ExtendAdd([[1,],[2,]], [[3,],[4,]]) = [[4,],[6,]]
    # which is different from np.concatenate()

    # record the indices that have different dimensions in array_A and array_B
    list_index_diff_dim, list_diff_dim_A, list_diff_dim_B = GetIndexDiffDim(array_A, array_B)
    # extend array_A with 0 on the indices that have different dimensions in array_B
    # the number of extended dimension = the dimension of the same index of array_B
    array_A_extend0 = Extend0(array_A,
                              list_extend_index=list_index_diff_dim,
                              list_extend_dim=list_diff_dim_B,
                              extend_to="tail",
                              )
    # extend array_B with 0 on the indices that have different dimensions in array_A
    # the number of extended dimension = the dimension of the same index of array_A
    array_B_extend0 = Extend0(array_B,
                              list_extend_index=list_index_diff_dim,
                              list_extend_dim=list_diff_dim_A,
                              extend_to="head",
                              )

    # check the shape
    if (array_A_extend0.shape != array_B_extend0.shape):
        print("----Error from ExtendAdd(): array_A_extend0.shape == array_B_extend0.shape violated !----")
        print(array_A_extend0.shape)
        print(array_B_extend0.shape)
        exit(1)
    else:
        # after extension, these two arrays have the same shape
        # which means that they can be directly added together
        array_C = array_A_extend0 + array_B_extend0
        return array_C

def GetIndexDiffDim(array_A, array_B):
    '''
    :param array_A: one of the two arrays to be compared.
    :param array_B: one of the two arrays to be compared.
    :return: a list of indices that have different dimensions in array_A and array_B,
    as well as two lists of these different dimensions for array_A and array_B respectively.
    '''
    # the two arrays are demanded to have the same order
    array_A = np.asarray(array_A)
    array_B = np.asarray(array_B)
    if (len(array_A.shape) != len(array_B.shape)):
        print("----Error from GetIndexDiffDim(): order_A == order_B violated !----")
        exit(1)
    # create lists to record the indices and corresponding dimensions
    list_index = []
    list_dim_A = []
    list_dim_B = []
    for index in range(0, len(array_A.shape)):
        dim_A_index = array_A.shape[index]
        dim_B_index = array_B.shape[index]
        # if this index has different dimensions for array_A and array_B, record it
        if (dim_A_index != dim_B_index):
            list_index.append(index)
            list_dim_A.append(dim_A_index)
            list_dim_B.append(dim_B_index)
    return list_index, list_dim_A, list_dim_B

def Extend0(array_A, list_extend_index, list_extend_dim, extend_to="tail"):
    '''
    :param array_A: array to extend.
    :param list_extend_index: list of the indices to extend.
    :param list_extend_dim: list of additional dimensions to extend.
    :param extend_to: {"large","small"} determine the extension direction.
    :return: array with extended zero elements.
    '''
    # judge legality
    if (len(list_extend_dim) != len(list_extend_index)):
        print("----Error from Extend0(): len(list_extend_dim) == len(list_extend_index) violated !----")
        exit(1)

    # extend zero elements = concatenate a zero-element array
    array_A = np.asarray(array_A)
    for i in range(0, len(list_extend_index)):
        extend_index = list_extend_index[i]
        extend_dim = list_extend_dim[i]

        # define a zero-element array to concatenate
        shape_0 = list(array_A.shape)
        shape_0[extend_index] = extend_dim
        array_0 = np.zeros(shape_0)

        # concatenate the given array with the zero-element array
        if (extend_to == "tail"):
            # "tail" means [1,2,3] ---> [1,2,3,0,]
            array_A = np.concatenate((array_A, array_0),
                                     axis=extend_index)
        elif (extend_to == "head"):
            # "head" means [1,2,3] ---> [0,1,2,3]
            array_A = np.concatenate((array_0, array_A),
                                     axis=extend_index)
        else:
            print("----Error: extend_to = head or tail !----")
            exit(1)

    return array_A

def CombineIndex(array_A, list_index_to_combine):
    '''
    :param array_A: array to combine.
    :param list_index_to_combine: list of combine indices.
    :return: an array obtained by combining the given indices of the original array.
    '''
    # process: A --transpose--> B --reshape--> C

    # A --transpose--> B
    # transpose in order to put the combine indices together
    transpose_index_A2B = sum(list_index_to_combine, [])
    # sum(list, [[1,2],[3,4]]) return [1,2,3,4]
    array_B = np.transpose(array_A, transpose_index_A2B)

    # get the array shape after combine from list_index_to_combine, i.e.'shape_C'
    shape_C = []
    for A_indices in list_index_to_combine:
        dim_of_C_index = 1
        for A_index in A_indices:
            dim_of_C_index = dim_of_C_index * array_A.shape[A_index]
        shape_C += [dim_of_C_index]

    # B --reshape--> C
    array_C = np.reshape(array_B, shape_C)

    return array_C

if __name__ == "__main__":

    # get the model parameters
    N, J, h, D, sweeps, cvg = GetParams()

    # get the Hamiltonian mpo (the variation mpo for the ground state)
    list_mpo_H = GetListHmpo(N, GetSiteHmpo_TFI(J, h))
    # get the ground state mps
    mps0 = mps(N, D, list_mpo_H)
    E0, counter = mps0.vMPS(sweeps=sweeps, cvg=cvg, if_print=True)

    # get the variation mpo for the first excited state
    list_mpo_H_excited = GetExcitedHmpo(list_mpo_H, mps0, E0)
    # get the first excited state mps
    mps1 = mps(N, D, list_mpo_H_excited)
    E1, counter = mps1.vMPS(sweeps=sweeps, cvg=cvg, if_print=True)

    # print the results compactly
    PrintE(E1, counter)

    # measure <\sigma_z>
    list_ave_sigma_z = mps1.Measure1SiteOpr(2 * GenSpinOpr('Sz'), range(0, N))
    print('<\sigma_z> =', list_ave_sigma_z)

    # measure <\sigma_x>
    list_ave_sigma_x = mps1.Measure1SiteOpr(2 * GenSpinOpr('Sx'), range(0, N))
    print('<\sigma_x> =', list_ave_sigma_x)

    # plot <\sigma_z> and <\sigma_x>
    PlotXYn(range(0, N),
            [list_ave_sigma_z,
             list_ave_sigma_x,
             ],
            list_ylabel=[LatexAve('\sigma^z'),
                         LatexAve('\sigma^x'),
                         ],
            name_x='x',
            name_y=LatexAve('\sigma'),
            )