#! /usr/bin/env python

# Copyright (c) 2015 William Lees

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

__author__ = 'Martin Rosellen'
__docformat__ = "restructuredtext en"

import argparse
import sys

import cairo
import pandas as pd

WIDTH, HEIGHT = 1500, 3000

RES_RADIUS = 60
COL_SPACING = 400
RES_Y_SPACING = 80
MARGIN = 100
FONT_SIZE = 24
DASH_SIZE = 10


def read_control_file(file_name):
    data_frame = pd.read_csv(file_name)
    n_chains = len(data_frame.Chain.unique())
    residues_control_file = data_frame.Id.unique()
    return data_frame, n_chains, residues_control_file


def read_decomp_table_file(file_name):
    data_frame = pd.read_csv(file_name, index_col=0)

    # remove rows not containing any energy values
    data_frame.dropna(axis=0, how='all', inplace=True)

    # omit first cell which contains 'Res'
    residues_decomp_table = data_frame.index

    return data_frame, residues_decomp_table


def read_hbonds_file(file_name, threshold):
    data_frame = pd.read_csv(file_name, names=['res1', 'res2', 'n_frames'])
    data_frame = data_frame[data_frame.n_frames > threshold]
    return data_frame


def check_residue_naming(control_file_residues, decomp_table_residues):
    control_residues_set = set(control_file_residues)
    decomp_residues_set = set(decomp_table_residues)

    common_residues = control_residues_set & decomp_residues_set

    if len(common_residues) == len(decomp_table_residues):
        print('Control file contains all elements of decomp table.')
    else:
        print('Missing residues in control file {}'.format(decomp_residues_set - control_residues_set))


def generate_residue_names_to_plot(control_file_data_frame, residues_decomp_table):
    selected_rows = control_file_data_frame[control_file_data_frame.Id.isin(residues_decomp_table)]

    residue_names = []
    for _, row in selected_rows.iterrows():
        residue_names.append('{}:{} {}'.format(row.Chain, row.Id.split()[0], row.Legend))
    return selected_rows, residue_names


def generate_column_x_coordinates(n_chains):
    middle = WIDTH / 2
    columns_x_coordinates = []

    if n_chains % 2 == 1:
        left_most_column_x = middle - ((n_chains // 2) * COL_SPACING)

        for c in range(n_chains):
            columns_x_coordinates.append(left_most_column_x + (c * COL_SPACING))
    else:
        left_most_column_x = middle - (COL_SPACING / 2) - (((n_chains / 2) - 1) * COL_SPACING)
        for c in range(n_chains):
            columns_x_coordinates.append(left_most_column_x + (c * COL_SPACING))

    return columns_x_coordinates


def main(args):
    parser = argparse.ArgumentParser(description='Plot residue-wise interaction energies.')
    parser.add_argument('control', help='The control file that determines which residues are plotted and how.')
    parser.add_argument('decomp', help='The decomp table produced by PairwiseDecompTable.')
    parser.add_argument('hbonds', help='The consolidated hbond file produced by ConsolidateHbonds.')
    parser.add_argument('thresh', help='The minimum threshold for hbonds (in how many frames should the hbond be '
                                       'present to be considered here).', default=150)
    parser.add_argument('output', help='File name of the diagram (PDF format).')
    parser.add_argument('summary', help='File name of the residue-wise energies table (CSV format).')
    parser.add_argument('-t', '--compare_thresh', help='Energy values must be higher than this threshold to be '
                                                       'considered (default 0.5 kcal/mol).', default=0.5)
    args = parser.parse_args()

    control_file_data_frame, n_chains, residues_control_file = read_control_file(args.control)

    decomp_table_data_frame, residues_decomp_table = read_decomp_table_file(args.decomp)

    check_residue_naming(residues_control_file, residues_decomp_table)

    hbonds_data_frame = read_hbonds_file(args.hbonds, int(args.thresh))

    # create cairo surface for plotting
    surface = cairo.PDFSurface(args.output, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    ctx.set_font_size(FONT_SIZE)

    selected_rows, residue_names_to_plot = generate_residue_names_to_plot(control_file_data_frame,
                                                                          residues_decomp_table)
    # calculate where the columns should be placed
    columns_x_coordinates = generate_column_x_coordinates(n_chains)


if __name__ == '__main__':
    main(sys.argv)
