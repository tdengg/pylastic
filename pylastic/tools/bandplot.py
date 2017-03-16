#!/usr/bin/env python

# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import sys
import numpy as np

try:
    import yaml
except ImportError:
    print("You need to install python-yaml.")
    sys.exit(1)
    
try:
    from yaml import CLoader as Loader
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from phonopy.units import VaspToTHz

def read_band_yaml(filename):
    data = yaml.load(open(filename), Loader=Loader)
    frequencies = []
    distances = []
    labels = []
    for j, v in enumerate(data['phonon']):
        if 'label' in v:
            labels.append(v['label'])
        else:
            labels.append(None)
        frequencies.append([f['frequency'] for f in v['band']])
        distances.append(v['distance'])

    return (np.array(distances),
            np.array(frequencies),
            data['segment_nqpoint'],
            labels)

def read_dos_dat(filename):
    dos = []
    frequencies = []
    for line in open(filename):
        if line.strip()[0] == '#':
            continue
        ary = [float(x) for x in line.split()]
        frequencies.append(ary.pop(0))
        dos.append(ary)
    return np.array(frequencies), np.array(dos)

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(band_labels=None,
                    dos=None,
                    dos_max=None, 
                    dos_min=None,
                    factor=1.0,
                    f_max=None, 
                    f_min=None,
                    is_gnuplot=False,
                    is_points=False,
                    is_vertial_line=False,
                    output_filename=None,
                    xlabel=None,
                    ylabel=None,
                    show_legend=False,
                    title=None)
parser.add_option("--band_labels", dest="band_labels", action="store", type="string",
                  help="Show labels at band segments")
parser.add_option("--dos", dest="dos", type="string",
                  help="Read dos.dat type file and plot with band structure")
parser.add_option("--dmax", dest="dos_max", type="float",
                  help="Maximum DOS plotted")
parser.add_option("--dmin", dest="dos_min", type="float",
                  help="Minimum DOS plotted")
parser.add_option("--factor", dest="factor", type="float",
                  help="Conversion factor to favorite frequency unit")
parser.add_option("--fmax", dest="f_max", type="float",
                  help="Maximum frequency plotted")
parser.add_option("--fmin", dest="f_min", type="float",
                  help="Minimum frequency plotted")
parser.add_option("--gnuplot", dest="is_gnuplot", action="store_true",
                  help="Output in gnuplot data style")
parser.add_option("--legend", dest="show_legend",
                  action="store_true",
                  help="Show legend")
parser.add_option("--line", "-l", dest="is_vertial_line",
                  action="store_true",
                  help="Vertial line is drawn at between paths")
parser.add_option("-o", "--output", dest="output_filename",
                  action="store", type="string",
                  help="Output filename of PDF plot")
parser.add_option("--xlabel", dest="xlabel", action="store", type="string",
                  help="Specify x-label")
parser.add_option("--ylabel", dest="ylabel", action="store", type="string",
                  help="Specify y-label")
parser.add_option("--points", dest="is_points",
                  action="store_true",
                  help="Draw points")
parser.add_option("-t", "--title", dest="title", action="store",
                  type="string", help="Title of plot")
(options, args) = parser.parse_args()

if options.output_filename:
    import matplotlib
    matplotlib.use('Agg')            

if not options.is_gnuplot:
    import matplotlib.pyplot as plt
    if options.band_labels:
        from matplotlib import rc
        rc('text', usetex=True)
    if options.dos:
        import matplotlib.gridspec as gridspec
        plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
        ax1 = plt.subplot(gs[0, 0])

colors = ['r-', 'b--', 'g-', 'c--', 'm-', 'k--', 'y-', 
          'r--', 'b--', 'g--', 'c--', 'm--', 'y--', 'k--']
if options.is_points:
    colors = [x + 'o' for x in colors]

count = 0


if len(args) == 0:
    filenames = ['band.yaml']
else:
    filenames = args

if options.is_gnuplot:
    print("# distance  frequency (bands are separated by blank lines)")

if options.dos:
    dos_frequencies, dos = read_dos_dat(options.dos)

for i, filename in enumerate(filenames):
    (distances,
     frequencies,
     segment_nqpoint,
     labels) = read_band_yaml(filename)

    end_points = [0,]
    for nq in segment_nqpoint:
        end_points.append(nq + end_points[-1])
    end_points[-1] -= 1
    segment_positions = distances[end_points]

    if all(x is None for x in labels):
        labels_at_ends = None
    else:
        labels_at_ends = [r"$%s$"%labels[n] for n in end_points]

    if options.is_gnuplot:
        print("# End points of segments: ")
        print("#   " + "%10.8f " * len(segment_positions) %
              tuple(segment_positions))
    elif options.is_vertial_line and len(filenames) == 1:
        for v in segment_positions[1:-1]:
            plt.axvline(x=v, linewidth=0.5, color='b')

    if options.is_gnuplot:
        for j, freqs in enumerate(frequencies.T):
            q = 0
            for nq in segment_nqpoint:
                for d, f in zip(distances[q:(q + nq)],
                                freqs[q:(q + nq)] * options.factor):
                    print("%f %f" % (d, f))
                q += nq
                print('')
            print('')
    else:
        q = 0
        for j, nq in enumerate(segment_nqpoint):
            if j == 0:
                plt.plot(distances[q:(q + nq)],
                         frequencies[q:(q + nq)] * options.factor,
                         colors[i])
            else:
                plt.plot(distances[q:(q + nq)],
                         frequencies[q:(q + nq)] * options.factor,
                         colors[i])
            q += nq

    if options.is_gnuplot:
        print('')
    
    plt.plot([],[],colors[i], label=filename)
    
if not options.is_gnuplot:
    if options.xlabel is None:
        plt.xlabel('Wave vector')
    else:
        plt.xlabel(options.xlabel)
    if options.ylabel is None:
        plt.ylabel('Frequency  (THz)')
    else:
        plt.ylabel(options.ylabel)
    
    plt.xlim(distances[0], distances[-1])
    if options.f_max is not None:
        plt.ylim(ymax=options.f_max)
    if options.f_min is not None:
        plt.ylim(ymin=options.f_min)
    else:
        plt.ylim(ymin=0)
    plt.axhline(y=0, linestyle=':', linewidth=0.5, color='b')
    if len(filenames) < 10:
        xticks = segment_positions
        if options.band_labels:
            band_labels = [x for x in options.band_labels.split()]
            if len(band_labels) == len(xticks):
                plt.xticks(xticks, band_labels)
            else:
                print("Numbers of labels and band segments don't match.")
                sys.exit(1)
        elif labels_at_ends:
            plt.xticks(xticks, labels_at_ends)
        else:
            plt.xticks(xticks, [''] * len(xticks))
    else:
        plt.xticks([])
    
    if options.title is not None:
        plt.title(options.title)

    if options.show_legend:
        plt.legend()

    if options.dos:
        ax2 = plt.subplot(gs[0, 1], sharey=ax1)
        plt.subplots_adjust(wspace=0.03)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.plot(dos, dos_frequencies)
        plt.xlabel('DOS')
        if options.f_max is not None:
            plt.ylim(ymax=options.f_max)
        if options.f_min is not None:
            plt.ylim(ymin=options.f_min)
        else:
            plt.ylim(ymin=0)
        if options.dos_max is not None:
            plt.xlim(xmax=options.dos_max)
        if options.dos_min is not None:
            plt.xlim(xmin=options.dos_min)
    
    if options.output_filename is not None:
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.savefig(options.output_filename)
    else:
        plt.show()
