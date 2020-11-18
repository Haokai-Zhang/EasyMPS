#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <plot_xy.py> contains functions of 2d-plot in a certain style.
'''

from matplotlib import pyplot as plt

def PlotXY(list_x, list_y, name_x='x', name_y='y', name_title=''):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y: list of scalars as data of y-axis
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :return: the graph of (x, y)
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    plt.rcParams['font.sans-serif'] = ['Arial']

    # label font and size
    font_xy = {'family': 'Arial', 'size': 19}
    font_title = {'family': 'Arial', 'size': 17}
    plt.xlabel(name_x, font_xy)
    plt.ylabel(name_y, font_xy)
    plt.title(name_title, font_title)

    # plot data
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)
    plt.plot(list_x,
             list_y,
             '-s',
             color=color_brick_red,
             markersize=9,
             markerfacecolor='none',
             markeredgewidth=1.6,
             linewidth=1.3,
             )

    # axis range and font
    number_font_size = 17
    plt.xticks(fontproperties='Arial', size=number_font_size)
    plt.yticks(fontproperties='Arial', size=number_font_size)

    # set space between picture and display box
    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.14, top=0.96)

    # border line width
    ax = plt.gca()
    axis_line_width = 1
    ax.spines['top'].set_linewidth(axis_line_width)
    ax.spines['bottom'].set_linewidth(axis_line_width)
    ax.spines['left'].set_linewidth(axis_line_width)
    ax.spines['right'].set_linewidth(axis_line_width)
    ax.tick_params(axis='x', direction='in', width=axis_line_width, length=5 * axis_line_width)
    ax.tick_params(axis='y', direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    # show the figure
    plt.show()

def PlotXYn(list_x, list_yn, list_ylabel=None,
            name_x='x', name_y='y', name_title='',
            list_color=None, list_marker=None,
            marker_size=9, marker_edge_width=1.6, line_width=1.3,
            x_lim=None, y_lim=None,
            ):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_yn: list of lists of scalars as data of y-axis
    :param list_ylabel: list of labels of different sets of data
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param list_color: list of marker color
    :param list_marker: list of marker shape
    :param marker_size: marker size
    :param marker_edge_width: marker edge width
    :param line_width: line width
    :param x_lim: x-axis limit
    :param y_lim: y-axis limit
    :return: the graph of (x, y1), (x, y2), ... 
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    plt.rcParams['font.sans-serif'] = ['Arial']

    # define preferred colors and markers
    if (list_color == None):
        list_color = [
            (204 / 255, 71 / 255, 38 / 255),
            (55 / 255, 110 / 255, 184 / 255),
            (222 / 255, 158 / 255, 42 / 255),
            ]
    if (list_marker == None):
        list_marker = [
            '-s',
            '-o',
            '-^',
        ]

    # label font and size
    font_title = {'family': 'Arial', 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': 'Arial', 'size': 19}
    plt.xlabel(name_x, font_xy, labelpad=-1)
    plt.ylabel(name_y, font_xy, labelpad=-1)

    # plot data
    list_handle = []
    for index_y in range(0, len(list_yn)):
        handle, = plt.plot(list_x,
                           list_yn[index_y],
                           list_marker[index_y % len(list_marker)],
                           color=list_color[index_y % len(list_color)],
                           markersize=marker_size,
                           markerfacecolor='none',
                           markeredgewidth=marker_edge_width,
                           linewidth=line_width,
                           )
        list_handle.append(handle)

    # set legend
    font_legend = {'family': 'Arial', 'size': 16}
    if (list_ylabel != None):
        plt.legend(handles=list_handle,
                   labels=list_ylabel,
                   loc='best',
                   prop=font_legend,
                   ncol=2,
                   frameon=False,
                   )

    # axis range and font
    tick_number_font_size = 17
    plt.xticks(fontproperties='Arial', size=tick_number_font_size)
    plt.yticks(fontproperties='Arial', size=tick_number_font_size)

    # set lim
    if (x_lim != None):
        plt.xlim(x_lim)
    if (y_lim != None):
        plt.ylim(y_lim)

    # set space between picture and display box
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.12, top=0.96)

    # border line width
    ax = plt.gca()
    axis_line_width = 1
    ax.spines['top'].set_linewidth(axis_line_width)
    ax.spines['bottom'].set_linewidth(axis_line_width)
    ax.spines['left'].set_linewidth(axis_line_width)
    ax.spines['right'].set_linewidth(axis_line_width)
    ax.tick_params(axis='x', direction='in', width=axis_line_width, length=5 * axis_line_width)
    ax.tick_params(axis='y', direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))

    # show the figure
    plt.show()

def PlotXnYn(list_xn, list_yn, list_ylabel=None,
            name_x='x', name_y='y', name_title='',
            list_color=None, list_marker=None,
            marker_size=9, marker_edge_width=1.6, line_width=1.3,
            x_lim=None, y_lim=None,
            ):
    '''
    :param list_xn: list of lists of scalars as data of x-axis
    :param list_yn: list of lists of scalars as data of y-axis
    :param list_ylabel: list of labels of different sets of data
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param list_color: list of marker color
    :param list_marker: list of marker shape
    :param marker_size: marker size
    :param marker_edge_width: marker edge width
    :param line_width: line width
    :param x_lim: x-axis limit
    :param y_lim: y-axis limit
    :return: the graph of (x, y1), (x, y2), ...
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    plt.rcParams['font.sans-serif'] = ['Arial']

    # define preferred colors and markers
    if (list_color == None):
        list_color = [
            (204 / 255, 71 / 255, 38 / 255),
            (55 / 255, 110 / 255, 184 / 255),
            (222 / 255, 158 / 255, 42 / 255),
            ]
    if (list_marker == None):
        list_marker = [
            '-s',
            '-o',
            '-^',
        ]

    # label font and size
    font_title = {'family': 'Arial', 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': 'Arial', 'size': 19}
    plt.xlabel(name_x, font_xy, labelpad=-1)
    plt.ylabel(name_y, font_xy, labelpad=-1)

    # plot data
    list_handle = []
    for index_y in range(0, len(list_yn)):
        handle, = plt.plot(list_xn[index_y],
                           list_yn[index_y],
                           list_marker[index_y % len(list_marker)],
                           color=list_color[index_y % len(list_color)],
                           markersize=marker_size,
                           markerfacecolor='none',
                           markeredgewidth=marker_edge_width,
                           linewidth=line_width,
                           )
        list_handle.append(handle)

    # set legend
    font_legend = {'family': 'Arial', 'size': 16}
    if (list_ylabel != None):
        plt.legend(handles=list_handle,
                   labels=list_ylabel,
                   loc='best',
                   prop=font_legend,
                   ncol=2,
                   frameon=False,
                   )

    # axis range and font
    tick_number_font_size = 17
    plt.xticks(fontproperties='Arial', size=tick_number_font_size)
    plt.yticks(fontproperties='Arial', size=tick_number_font_size)

    # set lim
    if (x_lim != None):
        plt.xlim(x_lim)
    if (y_lim != None):
        plt.ylim(y_lim)

    # set space between picture and display box
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.12, top=0.96)

    # border line width
    ax = plt.gca()
    axis_line_width = 1
    ax.spines['top'].set_linewidth(axis_line_width)
    ax.spines['bottom'].set_linewidth(axis_line_width)
    ax.spines['left'].set_linewidth(axis_line_width)
    ax.spines['right'].set_linewidth(axis_line_width)
    ax.tick_params(axis='x', direction='in', width=axis_line_width, length=5 * axis_line_width)
    ax.tick_params(axis='y', direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))

    # show the figure
    plt.show()

def PlotTwin(list_x, list_y1, list_y2, name_x='x', name_y1='y1', name_y2='y2', name_title=''):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y1: list of scalars as data of y1-axis
    :param list_y2: list of scalars as data of y1-axis
    :param name_x: name string of x-axis
    :param name_y1: name string of y1-axis
    :param name_y2: name string of y2-axis
    :param name_title: name string of figure title
    :return: the graph of (x, y1) and (x, y2) on twin-axis
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    plt.rcParams['font.sans-serif'] = ['Arial']

    # define preferred colors
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)

    # label font and size
    font_title = {'family': 'Arial', 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': 'Arial', 'size': 19}
    plt.xlabel(name_x, font_xy)
    color_y1 = color_brick_red
    plt.ylabel(name_y1, font_xy, color=color_y1)

    # plot y1 data
    plt.plot(list_x,
             list_y1,
             '-s',
             color=color_y1,
             markersize=9,
             markerfacecolor='none',
             markeredgewidth=1.6,
             linewidth=1.3,
             )

    # axis range and font
    tick_number_font_size = 17
    plt.xticks(fontproperties='Arial', size=tick_number_font_size)
    plt.yticks(fontproperties='Arial', size=tick_number_font_size, color=color_y1)

    # set space between picture and display box
    plt.subplots_adjust(left=0.12, right=0.88, bottom=0.14, top=0.96)

    # border line width
    ax = plt.gca()
    axis_line_width = 1
    ax.spines['top'].set_linewidth(axis_line_width)
    ax.spines['bottom'].set_linewidth(axis_line_width)
    ax.spines['left'].set_linewidth(axis_line_width)
    ax.spines['right'].set_linewidth(axis_line_width)
    ax.tick_params(axis='x', direction='in', width=axis_line_width, length=5 * axis_line_width)
    ax.tick_params(axis='y', direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    # twin axis
    ax_twin = ax.twinx()
    color_y2 = color_dark_blue
    ax_twin.plot(list_x,
                 list_y2,
                 '-o',
                 color=color_y2,
                 markersize=9,
                 markerfacecolor='none',
                 markeredgewidth=1.6,
                 linewidth=1.3,
                 )
    ax_twin.set_ylabel(name_y2, font_xy, color=color_y2)
    ax_twin.tick_params(axis='y', labelsize=tick_number_font_size, labelcolor=color_y2,
                        direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    # show the figure
    plt.show()

def LatexAve(opr_str):
    '''
    :param opr_str: string of operator
    :return: string of expectation value in LaTeX form <opr>
    '''
    return '$\langle \mathregular{' + opr_str + '} \\rangle$'

def LatexForm(opr_str):
    '''
    :param opr_str: string of operator
    :return: string of in LaTeX form
    '''
    return '$ \mathregular{' + opr_str + '} $'