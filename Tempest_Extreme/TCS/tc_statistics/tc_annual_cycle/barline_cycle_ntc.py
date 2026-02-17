#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Tue May  4 17:02:20 EDT 2021
if __name__ == '__main__':
    from misc.timer import Timer
    tt = Timer(f'start {__file__}')
import sys, os.path, os, glob, datetime
import xarray as xr, numpy as np, pandas as pd, matplotlib.pyplot as plt
#more imports
#
if __name__ == '__main__':
    tt.check('end import')
#
#start from here
from wyfig01a_barline_cycle_ntc_pTang import wyplot as plot_cycle
from wyfig01b_bar_cycle_ntc_pTang_sharpness import wyplot as plot_sharpness
from wyfig01c_bar_cycle_ntc_pTang_dfrac import wyplot as plot_dfrac
 
if __name__ == '__main__':
    from wyconfig import * #my plot settings
    figsize = 6,5
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2)
    axes = []

    ax = fig.add_subplot(gs[1,0])
    plot_sharpness(ax=ax)
    axes.append(ax)

    ax = fig.add_subplot(gs[1,1])
    plot_dfrac(ax=ax)
    axes.append(ax)

    ax = fig.add_subplot(gs[0,:])
    plot_cycle(ax=ax)
    axes.append(ax)

    for ax,tag in zip(axes, list('bca')):
        ax.text(0, 1, f'{tag}  ', ha='right', va='bottom', transform=ax.transAxes, fontweight='bold', fontsize='large')
    
    #savefig
    if len(sys.argv) > 1 and sys.argv[1] == 'savefig':
        #figname = __file__.replace('.py', f'.png')
        figname = __file__.replace('.py', f'.pdf')
        wysavefig(figname)
    tt.check(f'**Done**')
    plt.show()
    
