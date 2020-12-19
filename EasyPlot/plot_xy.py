#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-10-30 16:20
Description: EasyMPS project. <plot_xy.py> contains functions of 2d-plot in a certain style.
'''

from matplotlib import pyplot as plt
import numpy as np

def PlotXY(list_x, list_y,
           name_x='x', name_y='y', name_title='',
           font_size_name_xy=19, font_size_title=17,
           marker_color=None, marker_shape=None, marker_face_color='none',
           marker_size=9, marker_edge_width=1.6, line_width=1.3,
           font_size_tick_label=17,
           x_lim=None, y_lim=None,
           x_power_limit=(-1, 3), y_power_limit=(-1, 3),
           legend_y=None, where_legend='best', font_size_legend=16,
           adjust_left=0.12, adjust_right=0.95, adjust_bottom=0.12, adjust_top=0.96,
           dash_y=None,
           ):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y: list of scalars as data of y-axis
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param font_size_name_xy: font size of x-y axis name
    :param font_size_title: font size of title
    :param marker_color: marker color
    :param marker_shape: marker shape
    :param marker_face_color: marker face color
    :param marker_size: marker size
    :param marker_edge_width: marker edge width
    :param line_width: line width
    :param font_size_tick_label: font size of tick label
    :param x_lim: x-axis limit
    :param y_lim: y-axis limit
    :param x_power_limit: x-axis power limits
    :param y_power_limit: y-axis power limits
    :param legend_y: legend of y-data
    :param where_legend: position of legend
    :param font_size_legend: font size of legend
    :param adjust_left: left spacing between axis and frame
    :param adjust_right: right spacing between axis and frame
    :param adjust_bottom: bottom spacing between axis and frame
    :param adjust_top: top spacing between axis and frame
    :param dash_y: y-value of horizontal dashed line
    :return: the graph of (x, y)
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    fig = plt.figure(figsize=(figure_width, figure_height))  # size in inches
    ax = fig.add_subplot(1, 1, 1)

    # set all fonts
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # set ticks in
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # define preferred colors and markers
    if (marker_color == None):
        marker_color = (204 / 255, 71 / 255, 38 / 255)
        #(55 / 255, 110 / 255, 184 / 255)
        #(222 / 255, 158 / 255, 42 / 255)
    if (marker_shape == None):
        marker_shape = '-s'
        #'-o'
        #'-^'

    # plot data
    handle, = ax.plot(
        list_x,
        list_y,
        marker_shape,
        color=marker_color,
        markersize=marker_size,
        markerfacecolor=marker_face_color,
        markeredgewidth=marker_edge_width,
        linewidth=line_width,
    )

    # label font and size
    font_xy = {'family': font_family, 'size': font_size_name_xy}
    font_title = {'family': font_family, 'size': font_size_title}
    plt.xlabel(name_x, font_xy)
    plt.ylabel(name_y, font_xy)
    plt.title(name_title, font_title)

    # set legend
    font_legend = {'family': font_family, 'size': font_size_legend}
    if (legend_y != None):
        plt.legend(
            handles=(handle,),
            labels=(legend_y,),
            loc=where_legend,
            prop=font_legend,
            ncol=2,
            frameon=False,
        )

    if (dash_y != None):
        ax.plot(
            list_x,
            [dash_y for i in range(0, len(list_x))],
            linestyle=(0, (7, 6)),
            color='grey',
            linewidth=line_width,
        )

    # axis range and font
    number_font_size = font_size_tick_label
    plt.xticks(fontproperties=font_family, size=number_font_size)
    plt.yticks(fontproperties=font_family, size=number_font_size)

    # set lim
    if (x_lim != None):
        plt.xlim(x_lim)
    if (y_lim != None):
        plt.ylim(y_lim)

    # set space between picture and display box
    plt.subplots_adjust(left=adjust_left,
                        right=adjust_right,
                        bottom=adjust_bottom,
                        top=adjust_top,
                        )

    # border line width
    ax = plt.gca()
    axis_line_width = 1
    ax.spines['top'].set_linewidth(axis_line_width)
    ax.spines['bottom'].set_linewidth(axis_line_width)
    ax.spines['left'].set_linewidth(axis_line_width)
    ax.spines['right'].set_linewidth(axis_line_width)
    ax.tick_params(axis='x', direction='in', width=axis_line_width, length=5 * axis_line_width)
    ax.tick_params(axis='y', direction='in', width=axis_line_width, length=5 * axis_line_width, right='on')

    ax.xaxis.get_major_formatter().set_powerlimits(x_power_limit)
    ax.yaxis.get_major_formatter().set_powerlimits(y_power_limit)
    ax.yaxis.get_offset_text().set_fontfamily(font_family)

    # show the figure
    plt.show()

def PlotXYn(list_x, list_yn, list_legend_y=None,
            name_x='x', name_y='y', name_title='',
            list_color=None, list_marker=None, list_marker_face_color=None,
            marker_size=9, marker_edge_width=1.6, line_width=1.3,
            num_column=2,
            x_lim=None, y_lim=None,
            ):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_yn: list of lists of scalars as data of y-axis
    :param list_legend_y: list of labels of different sets of data
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param list_color: list of marker color
    :param list_marker: list of marker shape
    :param list_marker_face_color: list of marker face color
    :param marker_size: marker size
    :param marker_edge_width: marker edge width
    :param line_width: line width
    :param num_column: number of legend columns
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
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # set ticks in
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'


    # define preferred colors and markers
    if (list_color == None):
        list_color = [
            #'black',
            'royalblue',
            'darkorange',
            'tomato',
            'darkgreen',
            'saddlebrown',
        ]
    if (list_marker == None):
        list_marker = [
            '-s',
            '->',
            '-o',
            '-<',
            '-v'
        ]
    if (list_marker_face_color == None):
        list_marker_face_color = [
            #'black',
            'royalblue',
            'darkorange',
            'tomato',
            'darkgreen',
            'saddlebrown',
        ]

    # label font and size
    font_title = {'family': font_family, 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': font_family, 'size': 19}
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
                           markerfacecolor=list_marker_face_color[index_y % len(list_color)],
                           markeredgewidth=marker_edge_width,
                           linewidth=line_width,
                           )
        list_handle.append(handle)

    # set legend
    font_legend = {'family': font_family, 'size': 16}
    if (list_legend_y != None):
        plt.legend(handles=list_handle,
                   labels=list_legend_y,
                   loc='best',
                   prop=font_legend,
                   ncol=num_column,
                   frameon=False,
                   )

    # axis range and font
    tick_number_font_size = 17
    plt.xticks(fontproperties=font_family, size=tick_number_font_size)
    plt.yticks(fontproperties=font_family, size=tick_number_font_size)

    # set lim
    if (x_lim != None):
        plt.xlim(x_lim)
    if (y_lim != None):
        plt.ylim(y_lim)

    # set space between picture and display box
    plt.subplots_adjust(left=0.11, right=0.95, bottom=0.12, top=0.96)

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
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

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
    font_title = {'family': font_family, 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': font_family, 'size': 19}
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
    font_legend = {'family': font_family, 'size': 16}
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
    plt.xticks(fontproperties=font_family, size=tick_number_font_size)
    plt.yticks(fontproperties=font_family, size=tick_number_font_size)

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

def PlotTwin(
        list_x1, list_y1,
        list_x2=None, list_y2=None,
        name_x='x', name_y1='y1', name_y2='y2', name_title='',
):
    '''
    :param list_x1: list of scalars as data of x-axis for y1
    :param list_y1: list of scalars as data of y1-axis
    :param list_x2: list of scalars as data of x-axis for y2
    :param list_y2: list of scalars as data of y2-axis
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
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # define preferred colors
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)

    # label font and size
    font_title = {'family': font_family, 'size': 17}
    plt.title(name_title, font_title)

    font_xy = {'family': font_family, 'size': 19}
    plt.xlabel(name_x, font_xy)
    color_y1 = color_brick_red
    plt.ylabel(name_y1, font_xy, color=color_y1)

    # plot y1 data
    plt.plot(list_x1,
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
    plt.xticks(fontproperties=font_family, size=tick_number_font_size)
    plt.yticks(fontproperties=font_family, size=tick_number_font_size, color=color_y1)

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
    if (list_x2 == None) or (list_y2 == None):
        list_x2 = []
        list_y2 = []
        for i in range(0, len(list_x1)-1):
            list_x2.append((list_x1[i + 1] + list_x1[i]) * 0.5)
            list_y2.append((list_y1[i + 1] - list_y1[i]) / (list_x1[i + 1] - list_x1[i]))
    ax_twin.plot(list_x2,
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

def PlotXLogY(list_x, list_y, name_x='x', name_y='y', name_title=''):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y: list of scalars as data of y-axis
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :return: the graph of (x, log(y))
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # label font and size
    font_xy = {'family': font_family, 'size': 20}
    font_title = {'family': font_family, 'size': 17}
    plt.xlabel(name_x, font_xy)
    plt.ylabel(name_y, font_xy)
    plt.title(name_title, font_title)

    # plot data
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)
    plt.semilogy()
    plt.plot(list_x,
             list_y,
             '-s',
             color=color_brick_red,
             markersize=9,
             markerfacecolor='none',
             markeredgewidth=1.6,
             linewidth=1.3,
             )

    fit_coeff = np.polyfit(list_x, np.log(list_y), 1)
    fit_slope = float(fit_coeff[0])
    fit_func = np.poly1d(fit_coeff)
    fit_value = np.exp(fit_func(list_x))
    fit_handle, = plt.plot(list_x,
                           fit_value,
                           linestyle=(0, (9, 12)),
                           color=color_brick_red,
                           linewidth=1,
                           alpha=0.6,
                           )
    plt.legend(handles=(fit_handle,),
               labels=(name_y + '$ \propto e^{- ' + name_x + '/' + str('%.2f' % (-1 / fit_slope)) + '}$',),
               loc='best',
               prop={'family': font_family, 'size': 18},
               ncol=1,
               frameon=False,
               )


    # axis range and font
    number_font_size = 17
    plt.xticks(fontproperties=font_family, size=number_font_size)
    plt.yticks(fontproperties=font_family, size=number_font_size)

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

def PlotLogXY(list_x, list_y,
              name_x='x', name_y='y', name_title='',
              head_cut=0,
              legend_y=None,
              marker_color=None, marker_shape=None, marker_face_color='none',
              marker_size=9, marker_edge_width=1.6, line_width=1.3,
              x_lim=None, y_lim=None,
              ):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y: list of scalars as data of y-axis
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param head_cut: number of data points to cast off from head in fitting
    :param legend_y: legend of y-data
    :param marker_color: marker color
    :param marker_shape: marker shape
    :param marker_face_color: marker face color
    :param marker_size: marker size
    :param marker_edge_width: marker edge width
    :param line_width: line width
    :param x_lim: x-axis limit
    :param y_lim: y-axis limit
    :return: the graph of (log(x), y)
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # set ticks in
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # define preferred colors and markers
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)
    if (marker_color == None):
        marker_color = color_brick_red
    if (marker_shape == None):
        marker_shape = '-s'
        #'-o'
        #'-^'

    # label font and size
    font_xy = {'family': font_family, 'size': 20}
    font_title = {'family': font_family, 'size': 17}
    plt.xlabel(name_x, font_xy)
    plt.ylabel(name_y, font_xy)
    plt.title(name_title, font_title)

    # plot data
    plt.semilogx()
    plt.plot(list_x,
             list_y,
             marker_shape,
             color=marker_color,
             markersize=marker_size,
             markerfacecolor=marker_face_color,
             markeredgewidth=marker_edge_width,
             linewidth=line_width,
             )

    fit_coeff = np.polyfit(np.log(list_x[head_cut:]), list_y[head_cut:], 1)
    fit_slope = float(fit_coeff[0])
    fit_func = np.poly1d(fit_coeff)
    fit_value = fit_func(np.log(list_x))
    fit_handle, = plt.plot(list_x,
                           fit_value,
                           linestyle=(0, (9, 12)),
                           color=marker_color,
                           linewidth=1,
                           alpha=0.6,
                           )
    if (legend_y == None):
        legend_y = name_y + '$ \propto' + str('%.2f' % (fit_slope)) + '~ln($' + name_x + '$)$'
    plt.legend(handles=(fit_handle,),
               labels=(legend_y,),
               loc='best',
               prop={'family': font_family, 'size': 18},
               ncol=1,
               frameon=False,
               )


    # axis range and font
    number_font_size = 17
    plt.xticks(fontproperties=font_family, size=number_font_size)
    plt.yticks(fontproperties=font_family, size=number_font_size)

    # set lim
    if (x_lim != None):
        plt.xlim(x_lim)
    if (y_lim != None):
        plt.ylim(y_lim)

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

def PlotLogXLogY(list_x, list_y,
                 name_x='x', name_y='y', name_title='',
                 head_cut=0):
    '''
    :param list_x: list of scalars as data of x-axis
    :param list_y: list of scalars as data of y-axis
    :param name_x: name string of x-axis
    :param name_y: name string of y-axis
    :param name_title: name string of figure title
    :param head_cut: number of data points to cast off from head in fitting
    :return: the graph of (log(x), log(y))
    '''

    # define figure size
    centimeter_to_inch = 1 / 2.54
    figure_width = centimeter_to_inch * 21
    figure_height = figure_width * 0.66
    plt.figure(figsize=(figure_width, figure_height))  # size in inches

    # set all fonts
    font_family = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = font_family

    # label font and size
    font_xy = {'family': font_family, 'size': 20}
    font_title = {'family': font_family, 'size': 17}
    plt.xlabel(name_x, font_xy)
    plt.ylabel(name_y, font_xy)
    plt.title(name_title, font_title)

    # plot data
    color_brick_red = (204 / 255, 71 / 255, 38 / 255)
    color_dark_blue = (55 / 255, 110 / 255, 184 / 255)
    color_dull_yellow = (222 / 255, 158 / 255, 42 / 255)
    plt.loglog()
    plt.plot(list_x,
             list_y,
             '-s',
             color=color_brick_red,
             markersize=9,
             markerfacecolor='none',
             markeredgewidth=1.6,
             linewidth=1.3,
             )

    fit_coeff = np.polyfit(np.log(list_x[head_cut:]), np.log(list_y[head_cut:]), 1)
    fit_slope = float(fit_coeff[0])
    fit_func = np.poly1d(fit_coeff)
    fit_value = np.exp(fit_func(np.log(list_x)))
    fit_handle, = plt.plot(list_x,
                           fit_value,
                           linestyle=(0, (9, 12)),
                           color=color_brick_red,
                           linewidth=1,
                           alpha=0.6,
                           )
    plt.legend(handles=(fit_handle,),
               labels=(name_y + '$ \propto ' + name_x + '^{' + str('%.2f' % (fit_slope)) + '}$',),
               loc='best',
               prop={'family': font_family, 'size': 18},
               ncol=1,
               frameon=False,
               )


    # axis range and font
    number_font_size = 17
    plt.xticks(fontproperties=font_family, size=number_font_size)
    plt.yticks(fontproperties=font_family, size=number_font_size)

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

def LatexAve(opr_str):
    '''
    :param opr_str: string of operator
    :return: string of expectation value in LaTeX form <opr>
    '''
    return '$\langle ' + opr_str + ' \\rangle$'

def LatexForm(opr_str):
    '''
    :param opr_str: string of operator
    :return: string of in LaTeX form
    '''
    return '$ ' + opr_str + ' $'