#!/usr/bin/env python
# coding: utf-8

# To set up pyvista use Miniconda and pip e.g.
#
#  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
#  chmod +x Miniconda3-latest-MacOSX-x86_64.sh
#  ./Miniconda3-latest-MacOSX-x86_64.sh -b -p miniconda3
# source miniconda3/bin/activate
# conda create --name vtk_env python=3.8
# conda activate vtk_env
# conda install nodejs matplotlib
# pip install jupyter pyvista ipyvtk_simple
# 
# can be be run as a script e.g.
#
# python plot-field.py
#
# or in notebook - 
#
# jupyter notebook

import vtk
import pyvista as pv
from matplotlib import cm
import numpy as np

fcode='θ'
# fcode='state_debug'

mesh22 = pv.read('Coupler_UnitTest_atmosphere_long_num0004.pvtu')

xc=mesh22.get_array('xc')
yc=mesh22.get_array('yc')
zc=mesh22.get_array('zc')


R=np.sqrt(np.square(xc)+np.square(yc)+np.square(zc))

# Mask to surface and "top" half of sphere.
tol=1e-1 
sval=( np.abs(R[0]-R)>tol ) | ( zc>0 )

l=np.array(list(range(len( R )))); l[ np.abs(R[0]-R)<tol ];len(l);len( R[np.abs(R[0]-R)<tol] )
ptnan=l[sval]
ptplot=l[~sval]

mesh22[fcode][ptnan]=np.nan

msurf=mesh22.extract_points(ptplot,adjacent_cells=True,include_cells=True)

plotter = pv.Plotter()
# plotter.add_mesh(msurf,style='points',show_edges=True,scalars='θ',culling=True,cmap='jet',lighting=True,nan_opacity=0,point_size=10)
plotter.add_mesh(msurf,style='surface',show_edges=False,scalars=fcode,culling=True,cmap='jet',lighting=True,nan_opacity=0)
# plotter.add_mesh(msurf,style='points',scalars='θ',culling=True,cmap='jet',lighting=True,nan_opacity=0)
plotter.show()
