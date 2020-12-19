#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import linecache
import re

def ReadNumberFromText(text_path):
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
    return linecache.getline(text_path, line_num).strip()

def ReadNumberFromString(str):
    list_num_elem = re.findall('([-+]?\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?', str)
    list_num_float = []
    for num_elem in list_num_elem:
        # E = -1.1e-5  --->   num_elem[0]='-1.1', num_elem[1]='.1', num_elem[2]='e-5'
        list_num_float.append(float((num_elem[0] + num_elem[2])))
    return list_num_float