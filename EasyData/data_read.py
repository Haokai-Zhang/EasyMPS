#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-20 14:00
Description: EasyMPS project. <data_read.py> contains functions to read data from text file.
'''

import linecache
import re

def ReadNumberFromText(text_path):
    '''
    :param text_path: path of text file to read.
    :return: list of data for each line.
    '''
    # count the line number
    total_lines = 0
    for index_line, line in enumerate(open(text_path, 'r')):
        total_lines += 1
    # read number from each line
    list_num_line = []
    for line_num in range(0, total_lines):
        current_line = ReadStringFromText(text_path, line_num)
        list_num_float = ReadNumberFromString(current_line)
        list_num_line.append(list_num_float)

    return list_num_line

def ReadStringFromText(text_path, line_num):
    '''
    :param text_path: path of text file to read.
    :param line_num: line number.
    :return: the string of corresponding line.
    '''
    return linecache.getline(text_path, line_num).strip()

def ReadNumberFromString(str):
    '''
    :param str: string to read.
    :return: list of number from this string.
    '''
    list_num_elem = re.findall('([-+]?\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?', str)
    list_num_float = []
    for num_elem in list_num_elem:
        # E = -1.1e-5  --->   num_elem[0]='-1.1', num_elem[1]='.1', num_elem[2]='e-5'
        list_num_float.append(float((num_elem[0] + num_elem[2])))
    return list_num_float