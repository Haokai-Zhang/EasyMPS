#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-07-04 11:30
Description: EasyMPS project. <class_named_data.py> contains a simple class of data with name.
'''

class named_data(object):

    def __init__(self, name_str, data):
        # give each data a name to quickly decide whether two of them are identical
        self.name = name_str
        self.data = data

    def __eq__(self, other):
        # equal if they have the same name
        return self.name == other.name