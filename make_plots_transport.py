'''
Run plots for drifters from these simulations.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
from glob import glob
import op
# from matplotlib.colors import LogNorm
# import prettyplotlib as ppl

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# Grid info
# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
loc = '/Users/kthyng/Documents/research/postdoc/'
grid = tracpy.inout.readgrid(loc, llcrnrlat=27.5, urcrnrlat=30.5, llcrnrlon=-97.6, urcrnrlon=-87.5)
# actually using psi grid here despite the name
xr = grid['xpsi']
yr = grid['ypsi']

# What combinations of plots do I want to do?
# * One plot for each month, with subplots the different years
# * One plot for the actual DWH timing, with subplots the different years
# * Average all years together and show months or quarters individually 

## ------ For each month, do subplots of years ----- ##

months = np.arange(1,13)
mname = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
            'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
years = np.arange(2004, 2011)

Smax = 0

for i, month in enumerate(months):
# month = months[1]

    S = []

    for j, year in enumerate(years):

        ## Read in transport info for all days that month/year combo ##
        Files = glob('tracks/' + str(year) + '-'
                    + str(month).zfill(2) + '-*gc.nc')

        U = 0; V = 0
        for File in Files:
            d = netCDF.Dataset(File)
            U += d.variables['U'][:]
            V += d.variables['V'][:]
            d.close

        # S is at cell centers, minus ghost points
        Stemp = np.sqrt(op.resize(U[:,1:-1],0)**2 + op.resize(V[1:-1,:],1)**2)
        S.append(Stemp)

        Smax = max((Smax,Stemp.max()))



    ## Plot ##
    fig = plt.figure(figsize=(12,9.5))#, dpi=150)
    for i in xrange(len(S)):

        ax = fig.add_subplot(4,2,i+1)
        ax.set_frame_on(False) # kind of like it without the box
        ax.text(0.85, 0.0, str(years[i]), transform=ax.transAxes, color='0.3')

        if i==0 or i==2 or i==4:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0])
        elif i==1 or i==3:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), merslabels=[0,0,0,0], parslabels=[0,0,0,0])
        elif i==5:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), parslabels=[0,0,0,0])
        elif i==6:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))

        C =  np.log(S[i]/Smax)
        C = np.ma.masked_where(np.isinf(C), C)
        C = np.ma.masked_where(np.isnan(C), C)
        mappable = ax.pcolormesh(xr, yr, C, cmap='Greys', vmin=-6, vmax=-3)

    # adjust subplots
    fig.subplots_adjust(left=0.05, bottom=0.03, right=0.98, top=0.99, wspace=0.06, hspace=0.01)

    # Colorbar in upper left corner
    cax = fig.add_axes([0.53, 0.15, 0.45, 0.015]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=12) 
    clim = cb.get_clim()
    cb.set_ticks(np.linspace(clim[0],clim[1],6))
    cb.set_ticklabels(np.exp(np.linspace(clim[0],clim[1],6)).round(3))
    cb.set_label('Surface transport from near DWH site in ' + mname[month-1], fontsize=16)

    fig.savefig('figures/transport/' + str(mname[month-1]) + 'lowres.jpg', dpi=50)#, bbox_inches='tight', dpi=100)
    plt.close()