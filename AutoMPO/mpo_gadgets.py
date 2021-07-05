#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2021-07-04 20:00
Description: EasyMPS project.
'''

def BoundSort(list_A, list_B):
    '''
    :param list_A: list to be sorted.
    :param list_B: list to be sorted.
    :return: two lists after sorting simultaneously by list_A from smallest to largest.
    '''
    list_A, list_B = (list(bound_element) for bound_element in zip(*sorted(zip(list_A, list_B))))
    # zip(): pack, zip(*): unpack
    return list_A, list_B

def IsOrdered(list_A):
    '''
    :param list_A: integer list to be checked
    :return: whether in order (small to large)
    '''
    for i in range(0, len(list_A) - 1):
        if (list_A[i] >= list_A[i + 1]):
            return False
    return True

def IsClose(a, b, abs_tol=1e-14):
    '''
    :param a: float a.
    :param b: float b.
    :param abs_tol: absolute tolerance
    :return: close or not
    '''
    return abs(a - b) < abs_tol