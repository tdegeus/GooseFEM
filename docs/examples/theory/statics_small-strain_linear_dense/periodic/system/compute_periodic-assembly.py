#!/usr/bin/env python3
"""compute_periodic-assembly.py
  Create system matrix or vector for document.

Usage:
  compute_periodic-assembly.py [options]

Options:
  -h --help           Show help.
  --original-matrix   Copy original matrix to clipboard.
  --periodic-matrix   Copy periodic matrix to clipboard.
  --original-vector   Copy original vector to clipboard.
  --periodic-vector   Copy periodic vector to clipboard.
"""

from docopt import docopt
import pyperclip
import numpy as np
import re
import copy

# ==================================================================================================

def convertVector(b):

  # number of columns
  N = len(b)

  # remove dangling starting "+", replace empty entries by "."
  body = '& '.join([re.sub(r'(\+)(.*)',r'\2',j)+' ' if len(j)>0 else '. ' for j in b])
  # add line with column numbers
  body = r'\hline'+'\n'+'& '.join([r'$%d$'%(i+1) for i in range(N)])+r'\\ \hline'+'\n'+body

  # add header and footer
  txt  = '\n'+r'\resizebox{287mm}{!}{%'
  txt += '\n'+r'\renewcommand{\arraystretch}{1.5}'
  txt += '\n'+r'\begin{tabular}{|'+'C{31mm}'*N+r'|}'+'\n'
  txt += body+r'\\ \hline'
  txt += '\n'+r'\end{tabular}'
  txt += '\n'+r'}'

  print('Text copied to clipboard')
  pyperclip.copy(txt)

  return txt

# ==================================================================================================

def convertMatrix(A):

  # number of columns
  N = len(A[0])

  # add column number, remove dangling starting "+", replace empty entries by "."
  body = ['& '.join(['$%d $'%(irow+1)]+[re.sub(r'(\+)(.*)',r'\2',j)+' ' if len(j)>0 else '. ' for j in i])+r'\\' for irow,i in enumerate(A)]
  body = '\n'.join(body)
  # add line with column numbers
  body = r'\hline'+'\n'+'&'+' & '.join([r'$%d$'%(i+1) for i in range(N)])+r'\\ \hline'+'\n'+body

  # add header and footer
  txt  = '\n'+r'\resizebox{287mm}{!}{%'
  txt += '\n'+r'\renewcommand{\arraystretch}{1.5}'
  txt += '\n'+r'\begin{tabular}{|c|'+'C{31mm}'*N+r'|}'+'\n'
  txt += body+r' \hline'
  txt += '\n'+r'\end{tabular}'
  txt += '\n'+r'}'

  print('Text copied to clipboard')
  pyperclip.copy(txt)

  return txt

# ==================================================================================================

# connectivity
# (change to list starting from zero, for computations)
conn = np.array([
  [ 1, 2, 6, 5],
  [ 2, 3, 7, 6],
  [ 3, 4, 8, 7],
  [ 5, 6,10, 9],
  [ 6, 7,11,10],
  [ 7, 8,12,11],
  [ 9,10,14,13],
  [10,11,15,14],
  [11,12,16,15],
])-1

# periodic node-pairs [independent,dependent]
# (change to list starting from zero, for computations)
periodic = np.array([
  [ 1, 4],
  [ 1,16],
  [ 1,13],
  [ 5, 8],
  [ 9,12],
  [ 3,15],
  [ 2,14],
])-1

# color of each element
col = ['blue','cyan','green','magenta','olive','purple','red','teal','violet']

# number of DOFs
ndof = np.max(conn)+1

# --------------------------------------------------------------------------------------------------

# renumber DOFs: [ independent , dependent ]
# - original numbering
i      = np.arange(ndof)
# - increase DOF-number of dependent DOFs with a sufficiently large number
i[periodic[:,1]] += ndof
# - towards new numbering: order of modified DOFs (the old order is retained within each group)
i      = np.argsort(i)
# - original numbering
dof    = np.arange(ndof)
# - new numbering: increasing numbers starting from zero
dof[i] = np.arange(ndof)
# - link dependent DOFs to their independent sister: dependent DOFs are assembled to independent ones
dof[periodic[:,1]] = dof[periodic[:,0]]
# - number of (in)dependent DOFs
nnd    = periodic.shape[0]
nni    = ndof-nnd

# --------------------------------------------------------------------------------------------------

# allocate empty "A" and "b"
A = [['' for j in range(ndof)] for i in range(ndof)]
b = [ ''                       for i in range(ndof)]

# regular assembly of "A" and "b":
# - loop over all elements
# - loop over all DOFs in each element
for ielem,elem in enumerate(conn):

  for i,idof in enumerate(elem):
    b[idof] += r'+\textcolor{%s}{$b_{%d}^{(%d)}$}'%(col[ielem],i+1,ielem+1)

  for i,idof in enumerate(elem):
    for j,jdof in enumerate(elem):
      A[idof][jdof] += r'+\textcolor{%s}{$A_{%d%d}^{(%d)}$}'%(col[ielem],i+1,j+1,ielem+1)

# --------------------------------------------------------------------------------------------------

# allocate periodic system
Aper = [['' for j in range(nni)] for i in range(nni)]
bper = [ ''                      for i in range(nni)]

# regular assembly of "A" and "b":
# - loop over all elements
# - loop over all DOFs in each element
for ielem,elem in enumerate(conn):

  for i,idof in enumerate(elem):
    bper[dof[idof]] += r'+\textcolor{%s}{$b_{%d}^{(%d)}$}'%(col[ielem],i+1,ielem+1)

  for i,idof in enumerate(elem):
    for j,jdof in enumerate(elem):
      Aper[dof[idof]][dof[jdof]] += r'+\textcolor{%s}{$A_{%d%d}^{(%d)}$}'%(col[ielem],i+1,j+1,ielem+1)

# ==================================================================================================

# parse command-line options/arguments
args = docopt(__doc__,version='0.0.1')

if args['--original-matrix']: convertMatrix(A)
if args['--periodic-matrix']: convertMatrix(Aper)
if args['--original-vector']: convertVector(b)
if args['--periodic-vector']: convertVector(bper)
