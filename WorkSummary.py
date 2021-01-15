# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib notebook
import glob
import jmkfigure

# %% [markdown]
# ## Writeup:
#
# - <file:///Users/jklymak/AbHillInterAnalysis/writeup/Notes.tex>
# - <drafts5://open?uuid=3E1289D2-AC38-45B5-9A7D-6863F178DBE6>

# %% [markdown]
# ## 3-km runs
#
# ### Isolated Bathymetry
#
# Here we make a run with just an isolated patch of bathymetry.  The issue here is that a very energetic Taylor cap is spun up.  It precesses around the patch as a relatively slow (I think) mode-Kelvin wave.  This dominates the energetics.  
#
#
#
# Issue is that a model run with an isolated patch of topography creates a Taylor cap.  This ends up dominating the 
#
# Try a spunup version that does not start impulsively.  Idea is to reduce the Talyor cap...
#
# Looks like the spinup was working, but doesn't particularly help?  
#

# %%
for td in ['reduceddata/Iso3kmlowU10Amp305f141B059Patch/twod.nc', 'reduceddata/Iso3kmlowU10Amp305f141B059PatchSu/twod.nc']:
    with xr.open_dataset(td) as ds:
        fig, ax = plt.subplots(3, 3, constrained_layout=True, sharex=True, sharey=True)
        for i in range(9):
            try:
                dd = ds.isel(time=i)
                axx=ax.flat[i]
                axx.pcolormesh(dd.XC/1e3, dd.YC/1e3, dd.ETAN[:-1,:-1]-dd.ETAN.mean(), vmin=-0.11, vmax=0.1, cmap='RdBu_r')
            except:
                pass


# %% [markdown]
#
#
# ### Runs with all bathy:
#
# Note that these have an issue that Eta becomes nonesense if we do things 
#
# Work as a function of time:

# %%
ds = glob.glob('reduceddata/Iso3k*Base*')
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, constrained_layout=True)
for nn, d in enumerate(ds):
    with xr.open_dataset(f'{d}/work.nc') as ds:
        pc = axs.flat[nn].pcolormesh(np.arange(37), ds.Z, ds.work.T[:-1, :-1], cmap='RdBu_r', shading='flat',
                                     vmin=-3.5e5, vmax=3.5e5, rasterized=True)
        axs.flat[nn].set_title(d[12:], fontsize=9)
        axs.flat[nn].set_ylim(-4000, 0)
    fig.colorbar(pc, ax=axs, shrink=0.6)
jmkfigure.jmkprint('BaseRunsDiss', 'WorkSummary.ipynb')

# %%
with xr.open_dataset('reduceddata/Iso3kmlowU10Amp305f141B000Base/twod.nc') as ds:
    #ds = ds.isel(XC=50)
    fig, ax = plt.subplots(3, 3, constrained_layout=True)
    for i in range(9):
        dd = ds.isel(time=i)
        axx = ax.flat[i]
        pc = axx.pcolormesh(dd.XG/1e3, dd.YG/1e3, dd.ETAN, cmap='RdBu_r')
        print(dd.ETAN.mean())
        fig.colorbar(pc, ax=axx)
        #dd.ETAN.plot(ax=ax.flat[i], cmap='RdBu_r')

# %%
with xr.open_dataset('reduceddata/Iso1kmlowU10Amp305f141B000Base/twod.nc') as ds:
    print(ds)
    #ds = ds.isel(XC=50)
    fig, ax = plt.subplots(3, 4)
    for i in range(9):
        dd = ds.isel(time=i)
        print(dd.time)
        dd.ETAN.plot(ax=ax.flat[i], cmap='RdBu_r')

# %%
print(ds)

# %%
