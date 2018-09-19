#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse

import numpy as np
#import matplotlib.pyplot as plt
import h5py

filedir = '/data/uclhc/uci/user/dantrim/TruthAnalysisR21/h5_files/'
ttbar_DR_file = '{}/truth_ntup_ttbar_Wt_DR.h5'.format(filedir)
ttbar_DS_file = '{}/truth_ntup_ttbar_Wt_DS.h5'.format(filedir)
#ttbar_file = '{}/truth_ntup_410009.h5'.format(filedir)
#wt_DR_file = '{}/truth_ntup_Wt_DR.h5'.format(filedir)
#wt_DS_file = '{}/truth_ntup_Wt_DS.h5'.format(filedir)
wwbb_file = '{}/truth_ntup_444444.h5'.format(filedir)

def import_pyplot(on_brick = True) :
    if not on_brick :
        import matplotlib.pyplot as plt
        return plt
    elif on_brick :
        import matplotlib
        matplotlib.use("pdf")
        import matplotlib.pyplot as plt
        return plt

plt = import_pyplot()
from matplotlib.lines import Line2D

class Sample :
    def __init__(self, name = '', filename = '', color = '') :
        self.name = name
        self.filename = filename
        self.color = color

def chunks(input_h5_dataset, chunksize = 100000) :
    for x in range(0, input_h5_dataset.size, chunksize) :
        yield input_h5_dataset[x:x+chunksize]

def valid_idx(array) :
    lo = array > -np.inf
    hi = array < np.inf
    return lo & hi

def variables() :

    var = {}
    var['NN_d_wt'] = [0.5, 0, 7]
    return var

def var_names() :

    names = {}
    names['NN_d_wt'] = '$d_{Wt}$'
    return names

def make_plots(samples, args) :

    var_dict = variables()
    if args.var != '' :
        if args.var not in var_dict :
            print('ERROR requested variable is not in list of variables')
            sys.exit()
        tmp = {}
        tmp[args.var] = var_dict[args.var]
        var_dict = tmp

    variables_to_plot = var_dict.keys()
    n_vars = len(variables_to_plot)

    histo_data = {}
    weight_data = {}

    for s in samples :

        histo_data[s.name] = {}
        weight_data[s.name] = []

        for variable in variables_to_plot :
            histo_data[s.name][variable] = []

        with h5py.File(s.filename, 'r', libver = 'latest') as sample_file :
            dataset = sample_file['truth']
            for chunk in chunks(dataset) :

                # select only EM or ME events
                sel_idx = (chunk['isEM'] == 1) | (chunk['isME'] == 1)
                chunk = chunk[sel_idx]

                for ivar, variable in enumerate(variables_to_plot) :

                    weights = chunk['eventweight']
                    data = chunk[variable]

                    histo_data[s.name][variable].extend(data)
                    if ivar == 0 :
                        weight_data[s.name].extend(weights)

    # make plots
    for ivar, variable in enumerate(variables_to_plot) :
        print('[{}/{}] plotting {}'.format(ivar+1,n_vars, variable))

        fig, ax = plt.subplots(1,1)
        ax.tick_params(axis = 'both', which = 'both', labelsize = 16, direction = "in",
                    labelleft=True, bottom = True, top = True, right = True, left = True)
        which_grid = 'both'
        ax.grid(color = 'k', which = which_grid, linestyle = '--', lw = 1, alpha = 0.1)

        bound_info = var_dict[variable]
        bw = bound_info[0]
        x_lo = bound_info[1]
        x_hi = bound_info[2]
        bin_edges = np.arange(x_lo, x_hi+bw, bw)

        ax.set_xlabel(var_names()[variable], horizontalalignment = 'right', x = 1.0)
        ax.set_ylabel('a.u.', horizontalalignment = 'right', y = 1.0)

        hist_bin_content = {}
        for sample in samples :

            h, _ = ax.hist(histo_data[sample.name][variable], weights = weight_data[sample.name], color = sample.color, bins = bin_edges, label = sample.name, normed = True, histtype = 'step')
            hist_bin_content[sample.name] = h

        ax.set_yscale('log')
        handles = [Line2D([0], [0], color=s.color) for s in samples]
        labels = [s.name for s in samples]
        ax.legend(handles, labels, loc = 'best', frameon = False)

        fig.savefig('test_tt_wt_{}.pdf'.format(variable), bbox_inches = 'tight', dpi = 200)

    


    



def main() :

    parser = argparse.ArgumentParser(description = 'A script to plot truth-level comparisons between various ttbar and Wt samples')
    parser.add_argument('-v', '--var', help = 'Request a specific variable to plot', default = '')
    args = parser.parse_args()

    #ttbar_sample = Sample('ttbar', ttbar_file, 'r')
    #wt_dr_sample = Sample('wt_DR', wt_DR_file, 'b')
    #wt_ds_sample = Sample('wt_DS', wt_DS_file, 'g')
    ttbar_DR_sample = Sample('ttbar_DR', ttbar_DR_file, 'r')
    ttbar_DS_sample = Sample('ttbar_DS', ttbar_DS_file, 'b')
    wwbb_sample = Sample('wwbb', wwbb_file, 'g')

    samples = [ttbar_DR_sample, ttbar_DS_sample, wwbb_sample]

    make_plots(samples, args)
    

    

if __name__ == "__main__" :
    main()
