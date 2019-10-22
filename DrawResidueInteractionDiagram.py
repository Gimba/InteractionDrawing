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
import math
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

res_codes = {}
res_codes['ALA'] = 'A'
res_codes['ARG'] = 'R'
res_codes['ASN'] = 'N'
res_codes['ASP'] = 'D'
res_codes['CYS'] = 'C'
res_codes['CYX'] = 'C'
res_codes['GLU'] = 'E'
res_codes['GLN'] = 'Q'
res_codes['GLY'] = 'G'
res_codes['HIS'] = 'H'
res_codes['HIE'] = 'H'
res_codes['HID'] = 'H'
res_codes['HIP'] = 'H'
res_codes['ILE'] = 'I'
res_codes['LEU'] = 'L'
res_codes['LYS'] = 'K'
res_codes['MET'] = 'M'
res_codes['PHE'] = 'F'
res_codes['PRO'] = 'P'
res_codes['SER'] = 'S'
res_codes['THR'] = 'T'
res_codes['TRP'] = 'W'
res_codes['TYR'] = 'Y'
res_codes['VAL'] = 'V'

# hydrophobic colors
res_colours = {}
res_colours['FIWLVM'] = [1, 0, 0]
res_colours['YCA'] = [1, 0.5, 0.5]
res_colours['THGSQ'] = [0, 1, 0]
res_colours['RKNEPD'] = [0, 1, 1]


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
        amino_acid = res_codes[row.Id.split()[0]]
        residue_names.append('{}:{}{}'.format(row.Chain, amino_acid, row.Legend))
    return selected_rows, residue_names


def generate_column_x_coordinates(n_chains):
    x_middle = WIDTH / 2
    columns_x_coordinates = generate_coordinate_spread(x_middle, n_chains, COL_SPACING)

    return columns_x_coordinates


def generate_coordinate_spread(middle, n, spacing):
    coordinate_offsets = []
    if n % 2 == 1:
        left_most_column_x = middle - ((n // 2) * spacing)

        for c in range(n):
            coordinate_offsets.append(left_most_column_x + (c * spacing))
    else:
        left_most_column_x = middle - (spacing / 2) - (((n / 2) - 1) * spacing)
        for c in range(n):
            coordinate_offsets.append(left_most_column_x + (c * spacing))
    return coordinate_offsets


def generate_intra_column_y_coordinates(residues_in_each_chain):
    col_max = max(residues_in_each_chain.values())
    y_middle = HEIGHT / 2
    col_y_coordinates = {}

    for k, v in residues_in_each_chain.items():
        col_y_coordinates[k] = generate_coordinate_spread(y_middle, v, RES_Y_SPACING)

    return col_y_coordinates


def generate_residue_plotting_coordinates(n_chains, selected_rows):
    # calculate where the columns should be placed
    columns_x_coordinates = generate_column_x_coordinates(n_chains)

    residues_in_each_column = dict(selected_rows.Col.value_counts())

    intra_column_y_coordinates = generate_intra_column_y_coordinates(residues_in_each_column)

    residue_plotting_coordinates = {}
    for k, v in intra_column_y_coordinates.items():
        residue_plotting_coordinates[k] = [columns_x_coordinates[k - 1], v]

    return residue_plotting_coordinates


def generate_amino_acid_colour(amino_acid_code):
    colour = [v for k, v in res_colours.items() if amino_acid_code in k][0]
    return colour


def draw_residue(ctx, x, y, text):
    colour = generate_amino_acid_colour(text[2])
    ctx.set_source_rgb(colour[0], colour[1], colour[2])
    ctx.save()
    ctx.translate(x, y)
    ctx.scale(RES_RADIUS, RES_RADIUS / 2.)
    ctx.arc(0., 0., 1., 0., 2 * math.pi)
    ctx.restore()
    ctx.fill()
    ctx.set_source_rgb(0, 0, 0)
    (x_bearing, y_bearing, add_width, height, x_advance, y_advance) = ctx.text_extents("I")
    (x_bearing, y_bearing, width, height, x_advance, y_advance) = ctx.text_extents(text)
    ctx.move_to(x - (width + add_width) / 2, y + height / 2)
    ctx.show_text(text)


def plot_residues(residue_plotting_coordinates, residue_names_to_plot, chain_column_id_mapping, ctx):
    for col_id, coords in residue_plotting_coordinates.items():
        residue_names = [r for r in residue_names_to_plot if r[0] == chain_column_id_mapping[col_id]]
        residue_names = sorted(residue_names, key=lambda x: int(x[3:]))
        x = coords[0]
        for y, name in zip(coords[1], residue_names):
            draw_residue(ctx, x, y, name)


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

    # calculate where the residues should be plotted for every column
    residue_plotting_coordinates = generate_residue_plotting_coordinates(n_chains, selected_rows)

    chain_column_id_mapping = selected_rows.drop_duplicates('Chain')[['Chain', 'Col']].set_index('Col').to_dict()[
        'Chain']

    plot_residues(residue_plotting_coordinates, residue_names_to_plot, chain_column_id_mapping, ctx)


if __name__ == '__main__':
    main(sys.argv)
