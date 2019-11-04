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
import numpy as np
import pandas as pd

WIDTH, HEIGHT = 1500, 3000

RES_RADIUS = 60
COL_SPACING = 400
RES_Y_SPACING = 60
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


def read_decomp_table_file(file_name, cmpr, exclusive=""):
    data_frame = pd.read_csv(file_name, index_col=0)

    # remove rows not containing any energy values
    if not cmpr:
        data_frame.dropna(axis=0, how='all', inplace=True)
        data_frame.dropna(axis=1, how='all', inplace=True)

    # only keep column that contains selected residue
    if exclusive:
        data_frame = data_frame.loc[:, 'TYR  22']
        # since only one column is selected data_frame will be of type Series. This causes problems later and
        # therefore it get transformed into a dataframe object again
        data_frame = data_frame.to_frame()
        data_frame = data_frame.drop(data_frame[data_frame['TYR  22'] == 0].index)
    residues_decomp_table = data_frame.index

    return data_frame, residues_decomp_table


def get_residue_interaction_tuples(decomp_table_data_frame):
    interactions = {}

    # only iterate over lower triangle of data frame matrix - Not used
    # temp_data = decomp_table_data_frame.mask(np.triu(np.ones(decomp_table_data_frame.shape, dtype=np.bool_)))

    for row in decomp_table_data_frame.iterrows():
        dict_val = row[1].dropna().to_dict()
        if dict_val:
            interactions[row[0]] = row[1].dropna().to_dict()

    temp = {}
    for k0, v_dict in interactions.items():
        for k1, v1 in v_dict.items():
            temp[k0 + " " + k1] = v1

    return temp


def get_highlight_residues(residue_interaction_tuples, residue_to_highlight):
    residue_selection = []
    for k in residue_interaction_tuples.keys():
        res1 = k[0:7].strip()
        res2 = k[7:].strip()

        if residue_to_highlight.iloc[0].Id in [res1, res2]:
            residue_selection.append(res1)
            residue_selection.append(res2)

    return list(set(residue_selection))

def read_hbonds_file(file_name, threshold):
    data_frame = pd.read_csv(file_name, names=['res1', 'res2', 'n_frames'])
    data_frame = data_frame[data_frame.n_frames > threshold]
    data_frame.res1 = data_frame.res1 + ' ' + data_frame.res2
    data_frame.drop('res2', axis=1, inplace=True)
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
        residue_names.append([row.Id, '{}:{}{}'.format(row.Chain, amino_acid, row.Legend)])
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


def generate_residue_plotting_coordinates(n_chains, selected_rows, exclusive=False):
    # calculate where the columns should be placed
    columns_x_coordinates = generate_column_x_coordinates(n_chains)

    residues_in_each_column = dict(selected_rows.Col.value_counts())

    intra_column_y_coordinates = generate_intra_column_y_coordinates(residues_in_each_column)

    residue_plotting_coordinates = {}
    for k, v in intra_column_y_coordinates.items():
        residue_plotting_coordinates[k] = [[columns_x_coordinates[k - 1], l] for l in v]

    return residue_plotting_coordinates


def generate_amino_acid_colour(amino_acid_code):
    colour = [v for k, v in res_colours.items() if amino_acid_code in k][0]
    return colour


def change_residue_ids(control_file_data_frame, data_frame):
    control_file_data_frame['residue_naming_new'] = control_file_data_frame.apply(lambda r: r.Chain + ':' + r.Legend,
                                                                                  axis=1)
    try:
        df_index = np.asarray(data_frame.index)
    except:
        df_index = np.asarray([])
    try:
        df_columns = np.asarray(data_frame.columns)
    except:
        df_columns = np.asarray([])

    for row in control_file_data_frame.iterrows():
        if df_index.size > 0:
            df_index[df_index == row[1].Id] = row[1].residue_naming_new
        if df_columns.size > 0:
            df_columns[df_columns == row[1].Id] = row[1].residue_naming_new
        data_frame.replace(row[1].Id, row[1].residue_naming_new, inplace=True)
    data_frame.set_index(df_index)
    data_frame.columns = df_columns
    return data_frame


def draw_residue(ctx, x, y, text, annotation=False, highlight=False):
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

    if highlight:
        ctx.select_font_face("Arial",
                             cairo.FontSlant.NORMAL,
                             cairo.FontWeight.BOLD)
    else:
        ctx.select_font_face("Arial",
                             cairo.FontSlant.NORMAL,
                             cairo.FontWeight.NORMAL)
    ctx.show_text(text)

    if annotation:
        ctx.move_to(x + (RES_RADIUS * 0.35), y - (RES_RADIUS * 0.35))
        ctx.set_font_size(FONT_SIZE * 0.75)
        ctx.show_text('{0:+.2f}'.format(annotation))
        ctx.set_font_size(FONT_SIZE)

def plot_residues(residue_plotting_coordinates, residue_names_to_plot, chain_column_id_mapping, ctx,
                  annotate_list=[], residue_selection=[]):
    residue_coordinates = {}
    highlight = False
    for col_id, coords in residue_plotting_coordinates.items():
        residue_names = [r for r in residue_names_to_plot if r[1][0] == chain_column_id_mapping[col_id]]
        residue_names = sorted(residue_names, key=lambda x: int(x[1][3:]))
        # x = coords[0]
        if annotate_list:
            for coord, name, annotation in zip(coords, residue_names, annotate_list):
                x = coord[0]
                y = coord[1]
                if name[0] in residue_selection:
                    highlight = True
                draw_residue(ctx, x, y, name[1], annotation, highlight)
                residue_coordinates[name[0]] = [x, y]
                highlight = False
        else:
            for coord, name in zip(coords, residue_names):
                x = coord[0]
                y = coord[1]
                draw_residue(ctx, x, y, name[1])
                residue_coordinates[name[0]] = [x, y]

    return residue_coordinates


def plot_interactions(residue_interaction_tuples, residue_coordinates, ctx, hbonds, residue_to_highlight=False):
    hbonds = list(hbonds.res1)

    # added this ugly list here to check which interaction have already been painted
    painted = []

    # maintain a list that holds the residues that are in contact with the one we want to highlight
    residues_to_highlight = []
    for r1, energy in residue_interaction_tuples.items():

        res1 = r1[0:7].strip()
        res2 = r1[7:].strip()

        # continue if the flipped residue identifiers are in the ugly list
        if res2 + ' ' + res1 in painted:
            continue
        # skip ions
        if '-' in res2 or '-' in res1:
            continue
        painted.append(r1)

        coord1 = residue_coordinates[res1]
        coord2 = residue_coordinates[res2]
        ctx.move_to(coord1[0], coord1[1])
        if r1 in hbonds:
            ctx.set_source_rgb(1, 0, 0)
        else:
            ctx.set_source_rgb(0, 0, 0)

        if energy > 0:
            ctx.set_dash([DASH_SIZE])
        else:
            ctx.set_dash({})

        ctx.set_line_width(abs(energy))
        ctx.line_to(coord2[0], coord2[1])
        ctx.stroke()


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
    parser.add_argument('-c', '--compare_file', help='only display interactions that differ from those in this file')
    parser.add_argument('-a', '--annotate', help='annotate residues with energy values')
    parser.add_argument('-m', '--highlight_residue', help='highlight interactions of this residue')
    parser.add_argument('-e', action="store_true",
                        help='specify that only contacts to the highlighted residue should be plotted')

    args = parser.parse_args()

    args.compare_thresh = float(args.compare_thresh)

    control_file_data_frame, n_chains, residues_control_file = read_control_file(args.control)

    # check if residue to highlight is present in control file
    if args.highlight_residue:
        residue_to_highlight = control_file_data_frame[control_file_data_frame.Legend.str.contains(
            args.highlight_residue)]
    else:
        residue_to_highlight = False

    if not residue_to_highlight.iloc[0].Id and args.e:
        print('Plot contacts to highlight residue exclusively selected but residue not specified. Please specify '
              'residue using -m. Exiting')
        exit()

    decomp_table_data_frame, residues_decomp_table = read_decomp_table_file(args.decomp, args.compare_file, args.e)

    if args.compare_file:
        compare_file_data_frame, residues_compare_file = read_decomp_table_file(args.compare_file, args.compare_file)

        # check if index of compare file data frame contains a mutation of a residue on base of the residue
        # numbering, this code is to only replace mutated residues (does not really work properly so far,
        # since the assumption here is that both indices are in the same order and have the same length)
        if not compare_file_data_frame.index.equals(decomp_table_data_frame.index):
            cmp_index = list(compare_file_data_frame)
            decomp_index = list(decomp_table_data_frame.index)

            # check if residue ids differ
            for k, [cmp_id, decomp_id] in enumerate(zip(cmp_index, decomp_index)):
                if cmp_id != decomp_id:
                    print('Differing amino acid code for residue', cmp_id.split()[1], ':', cmp_id.split()[0],
                          '(comp file)', decomp_id.split()[0], '(decomp) --> will plot', decomp_id.split()[0])
                    cmp_index[k] = decomp_id

            # replace index and column headers
            compare_file_data_frame.index = cmp_index
            compare_file_data_frame.columns = cmp_index

        # substracting data frames (with fill_value=0 'nan's get set to zero for substraction)
        decomp_table_data_frame = decomp_table_data_frame.subtract(compare_file_data_frame, fill_value=0)

        # apply threshold
        decomp_table_data_frame = decomp_table_data_frame.mask(abs(decomp_table_data_frame) <= args.compare_thresh)

        # remove colmuns and rows only containing nan
        decomp_table_data_frame.dropna(axis=0, how='all', inplace=True)
        decomp_table_data_frame.dropna(axis=1, how='all', inplace=True)

        # redefine residues to plot
        residues_decomp_table = decomp_table_data_frame.index

    check_residue_naming(residues_control_file, residues_decomp_table)

    hbonds_data_frame = read_hbonds_file(args.hbonds, int(args.thresh))

    # create cairo surface for plotting
    surface = cairo.PDFSurface(args.output, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    ctx.set_font_size(FONT_SIZE)

    residue_interaction_tuples = get_residue_interaction_tuples(decomp_table_data_frame)

    if residue_to_highlight is not None:
        residue_selection = get_highlight_residues(residue_interaction_tuples, residue_to_highlight)
    else:
        residue_selection = False

    selected_rows, residue_names_to_plot = generate_residue_names_to_plot(control_file_data_frame,
                                                                          residues_decomp_table)

    # calculate where the residues should be plotted for every column
    residue_plotting_coordinates = generate_residue_plotting_coordinates(n_chains, selected_rows, args.e)

    chain_column_id_mapping = selected_rows.drop_duplicates('Chain')[['Chain', 'Col']].set_index('Col').to_dict()[
        'Chain']

    residue_coordinates = plot_residues(residue_plotting_coordinates, residue_names_to_plot, chain_column_id_mapping,
                                        ctx)

    plot_interactions(residue_interaction_tuples, residue_coordinates, ctx, hbonds_data_frame)

    # replotting residues so that they overlay the interaction lines
    if args.annotate:
        gains_losses = list(decomp_table_data_frame.agg(lambda x: x.sum()))
    else:
        gains_losses = False

    plot_residues(residue_plotting_coordinates, residue_names_to_plot, chain_column_id_mapping,
                  ctx, gains_losses)


if __name__ == '__main__':
    main(sys.argv)
