#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-12-20 14:00
Description: EasyMPS project. <data_read.py> contains functions to read data from text file.
'''

import linecache
import re
import json

def ReadJson(json_file_name):
    '''
    :param json_file_name: the name of file to read.
    :return: data read.
    '''
    with open(json_file_name + '.json', 'r') as file_data:
        list_data_read = json.load(file_data)
        file_data.close()
    return list_data_read

def WriteJson(list_write, json_file_name):
    '''
    :param list_write: data to write.
    :param json_file_name: the name of file to write.
    :return: None.
    '''
    with open(json_file_name + '.json', 'w') as file_save:
        json_data = json.dumps(list_write, indent=4)
        file_save.write(json_data)
        file_save.close()

def ReadNumberFromFile(file_path, line_start=0, line_end=None):
    '''
    :param file_path: path of text file to read.
    :return: list of data for each line.
    '''
    # count the line number
    total_lines = CountLines(file_path)
    # read number from each line
    if (line_end == None):
        line_end = total_lines
    list_num_for_each_line = []
    for line_num in range(line_start, line_end):
        list_num_float = ReadNumberFromLine(file_path, line_num)
        list_num_for_each_line.append(list_num_float)

    return list_num_for_each_line

def ReadNumberFromLine(file_path, line_num):
    '''
    :param file_path: path of the file to read.
    :param line_num: the number of the line to read.
    :return: the float number in the wanted line.
    '''
    list_num_float = ReadNumberFromString(ReadStringFromFile(file_path, line_num))
    return list_num_float

def CountLines(file_path):
    '''
    :param file_path: path of the file to read.
    :return: total line number.
    '''
    total_lines = 0
    for index_line, line in enumerate(open(file_path, 'r')):
        total_lines += 1
    return total_lines

def ReadStringFromFile(file_path, line_num):
    '''
    :param file_path: path of the file to read.
    :param line_num: line number.
    :return: the string of corresponding line.
    '''
    return linecache.getline(file_path, line_num).strip()

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