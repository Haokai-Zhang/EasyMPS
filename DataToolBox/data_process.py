#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-20 14:00
Description: EasyMPS project. <data_process.py> contains functions as primary tools to process data.
'''

import numpy as np

def SvdTrunc(mat, D_trunc=None, if_normalize=False):
    '''
    :param mat: matrix to be singular-value decomposed.
    :param D_trunc: number of singular values to keep, rest to truncate.
    :param if_normalize: {True, False} determine whether to normalize.
    :return: singular-value decomposition results with truncation error.
    '''
    # svd using <numpy linear algebra>
    U, S, VT = np.linalg.svd(mat, full_matrices=False)
    # truncate and evaluate the truncation error
    trunc_err = 0
    if (D_trunc != None) and (len(S) > D_trunc):
        trunc_err = (S[D_trunc : ] ** 2).sum()
        U = U[:, 0 : D_trunc]
        S = S[0 : D_trunc]
        VT = VT[0 : D_trunc, :]
        # normalize
        if (if_normalize == True):
            S /= np.sqrt(np.sum(S ** 2))
    return U, S, VT, trunc_err

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

def Coor2Idx(x, y, Ny):
    return x * Ny + y

def Idx2Coor(site, Ny):
    x = site // Ny
    y = site % Ny
    return x, y

def Idx2x(site, Ny):
    x = site // Ny
    return x

def Idx2y(site, Ny):
    y = site % Ny
    return y

def GenBondSquare2dNN(Nx, Ny, boundary_condition='OO'):
    '''
    :param Nx: 2d square lattice length
    :param Ny: 2d square lattice size
    :param boundary_condition: boundary condition = {'OO', 'OP', 'PO', 'PP'}
    # 'O': open,  'P': periodic,  'OP': x-open, y-periodic
    :return: list of 2-dimensional nearest neighbor bonds on square lattice
    '''

    # example:
    #  1---3---5---
    #  |   |   |
    #  0---2---4---

    list_hop = []
    if (boundary_condition == 'OO'):
        # 'OO': x-open, y-open
        for x in range(0, Nx):
            for y in range(0, Ny):
                site_1 = Coor2Idx(x, y, Ny)
                x_next_x = (x + 1) % Nx
                y_next_x = y
                x_next_y = x
                y_next_y = (y + 1) % Ny
                site_2_x = Coor2Idx(x_next_x, y_next_x, Ny)
                site_2_y = Coor2Idx(x_next_y, y_next_y, Ny)
                if (x < Nx - 1):
                    list_hop.append([ site_1, site_2_x ])
                if (y < Ny - 1):
                    list_hop.append([ site_1, site_2_y ])
    elif (boundary_condition == 'OP'):
        # 'OP': x-open, y-periodic
        for x in range(0, Nx):
            for y in range(0, Ny):
                site_1 = Coor2Idx(x, y, Ny)
                x_next_x = (x + 1) % Nx
                y_next_x = y
                x_next_y = x
                y_next_y = (y + 1) % Ny
                site_2_x = Coor2Idx(x_next_x, y_next_x, Ny)
                site_2_y = Coor2Idx(x_next_y, y_next_y, Ny)
                if (x < Nx - 1):
                    list_hop.append([ site_1, site_2_x ])
                if (y < Ny - 1) or (Ny > 2):
                    list_hop.append([ site_1, site_2_y ])
    elif (boundary_condition == 'PO'):
        # 'PO': x-periodic, y-open
        for x in range(0, Nx):
            for y in range(0, Ny):
                site_1 = Coor2Idx(x, y, Ny)
                x_next_x = (x + 1) % Nx
                y_next_x = y
                x_next_y = x
                y_next_y = (y + 1) % Ny
                site_2_x = Coor2Idx(x_next_x, y_next_x, Ny)
                site_2_y = Coor2Idx(x_next_y, y_next_y, Ny)
                if (x < Nx - 1) or (Nx > 2):
                    list_hop.append([ site_1, site_2_x ])
                if (y < Ny - 1):
                    list_hop.append([ site_1, site_2_y ])
    elif (boundary_condition == 'PP'):
        # 'PO': x-periodic, y-open
        for x in range(0, Nx):
            for y in range(0, Ny):
                site_1 = Coor2Idx(x, y, Ny)
                x_next_x = (x + 1) % Nx
                y_next_x = y
                x_next_y = x
                y_next_y = (y + 1) % Ny
                site_2_x = Coor2Idx(x_next_x, y_next_x, Ny)
                site_2_y = Coor2Idx(x_next_y, y_next_y, Ny)
                if (x < Nx - 1) or (Nx > 2):
                    list_hop.append([ site_1, site_2_x ])
                if (y < Ny - 1) or (Ny > 2):
                    list_hop.append([ site_1, site_2_y ])
    else:
        print('-----Error: boundary_condition = OO/OP/PO/PP !-----')
        exit()
    return list_hop

def BoundSort(list_A, list_B):
    '''
    :param list_A: list to be sorted.
    :param list_B: list to be sorted.
    :return: two lists after sorting simultaneously by list_A from smallest to largest.
    '''
    list_A, list_B = (list(bound_element) for bound_element in zip(*sorted(zip(list_A, list_B))))
    # zip(): pack, zip(*): unpack
    return list_A, list_B

def IsCvg(list_data, cvg, check_len):
    '''
    :param list_E: list of data.
    :param cvg: convergence tolerance.
    :param check_len: convergence check length.
    :return: convergence is {True, False}.
    '''
    flag_cvg = True
    if (len(list_data) < check_len):
        return False
    else:
        for index_data in range(0, check_len - 1):
            if (abs(list_data[-1 - index_data] - list_data[-1 - (index_data + 1)]) > cvg):
                flag_cvg = False
        return flag_cvg

def IsClose(a, b, abs_tol=1e-14):
    '''
    :param a: float a.
    :param b: float b.
    :param abs_tol: absolute tolerance
    :return: close or not
    '''
    return abs(a - b) < abs_tol

def ShiftDeg(x, list_y, prec=1e-1, shift_unit=8e-2):
    '''
    :param x: parameter of degenerate values.
    :param list_y: list of possibly degenerate values.
    :param prec: degenerate precision.
    :param shift_unit: unit length of shift along x to identify degenerate values.
    :return: list of slightly shifted parameters used to plot degenerate values in a distinguishable way.
    '''
    # record the shifted x
    list_x_shifted = []
    # check degenerate values
    index_y = 0
    while (index_y < len(list_y)):
        # probe forward
        deg = 1
        while (index_y + deg < len(list_y)) \
                and (abs(list_y[index_y + deg] - list_y[index_y]) < prec) :
            deg += 1
        # shift x according to these degenerate values
        for index_deg in range(0, deg):
            list_x_shifted.append(x - 0.5 * (deg - 1) * shift_unit + index_deg * shift_unit)
        # read from the next different value in the next step
        index_y = index_y + deg

    return list_x_shifted