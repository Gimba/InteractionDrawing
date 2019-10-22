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

import pandas as pd


def read_control_file(file_name):
    data_frame = pd.read_csv(file_name)
    chains = data_frame.Chain.unique()
    residues_control_file = data_frame.Id.unique()
    return data_frame, chains, residues_control_file


def read_decomp_table_file(file_name):
    data_frame = pd.read_csv(file_name, index_col=0)

    # remove rows not containing any energy values
    data_frame.dropna(axis=0, how='all', inplace=True)

    # omit first cell which contains 'Res'
    residues_decomp_table = data_frame.columns[:]
    
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

    control_file_data_frame, chains, residues_control_file = read_control_file(args.control)

    decomp_table_data_frame, residues_decomp_table = read_decomp_table_file(args.decomp)

    check_residue_naming(residues_control_file, residues_decomp_table)

    hbonds_data_frame = read_hbonds_file(args.hbonds, int(args.thresh))


if __name__ == '__main__':
    main(sys.argv)
