#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:43:56 2022

@author: willdesnoo
"""

import numpy as np
import sys
from pprint import pprint
from ase import Atoms
from ase.io import write

def read_log_lines(log_name):
    """
    

    Parameters
    ----------
    log_name : string
        name of input file, fmt=lammps dump file

    Returns
    -------
    latdict : dictionary of list of strings
        lattice parameters for each image
        step number is key, 3x3 list of strings
        
    datadict : dictionary of list of strings
        geometric cooordinates of image
        step number is key
        id type element x y z fx fy fz
    """
    dbool=False
    lbool=False
    nstruc=0
    datadict={}
    latdict={}
    with open(log_name) as log_f:  
        log_lines=log_f.readlines()
        log_dat = []
        latparams=[]
        ts=False
        for i,l in enumerate(log_lines):
            if  "ITEM: ATOMS" in l:
                dbool=True
                lbool=False
                continue
            
            if dbool:
                if "ITEM: TIMESTEP" in l:
                    ts=True
                    dbool=False
                    datadict[nstruc] = log_dat
                    log_dat=[]
                    latparams=[]
                    continue
                spl=l.split(' ')
                latdict[nstruc] = latparams
                log_dat.append(spl[:-1])
                
            if ts:
                nstruc=l[:-1]
                ts=False
                
            if "ITEM: BOX" in l:
                lbool=True
                continue
            
            if lbool:
                latparams.append([x for x in l.split(' ')])
                
    return latdict, datadict

def get_neb_ordered(data, stepi, stepf, steps):
    result=[]
    for i in range(stepi,stepf+1,steps):
        natoms=len(data[str(i)])
        od=reorder(data[str(i)])
        neb_inp=[' '.join([x[0], x[3], x[4], x[5]]) + '\n' for x in od]
        neb_inp.insert(0, str(natoms) + '\n')
        result.append(neb_inp)
    return result
    
def get_fmtd_lp(latparam):
    xlo,xhi,xy=', '.join(latparam[0]).split(', ')
    ylo,yhi,xz=', '.join(latparam[1]).split(', ')
    zlo,zhi,yz=', '.join(latparam[2]).split(', ')
    xlo = xlo.replace('e+01','e+00') # Fix a bug in lammps dump xlo
    return xlo,xhi,ylo,yhi,zlo,zhi,xy[:-1],xz[:-1],yz[:-1]

def write_initial_coords(dat,latparam):
    natoms=len(dat)
    ordl=reorder(dat)
    olines=[]
    # Adding beginning lines of inital_coords
    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz=get_fmtd_lp(latparam)
    olines.append("initial_coords.in (Written by md2neb.py) \n")
    olines.append("\n")
    olines.append(f"{str(natoms)} atoms \n")
    olines.append("5 atom types \n")
    olines.append(f"{xlo} {xhi} xlo xhi \n")
    olines.append(f"{ylo} {yhi} ylo yhi \n")
    olines.append(f"{zlo} {zhi} zlo zhi \n")
    olines.append(f"{xy} {xz} {yz} xy xz yz \n")
    olines.append("\n")
    olines.append("\n")
    olines.append("Atoms \n")
    olines.append("\n")
      
    for d in ordl:
        olines.append(" ".join([d[0], d[1], "0.0", d[3], d[4], d[5],"\n"]))
    with open ("initial_coords.in","w") as f:
        for l in olines:
            f.write(l)
            
            
def reorder(li):
    return sorted(li,key=lambda x: x[1])
    

def column(matrix, i):
    return [row[i] for row in matrix]


def main(argv):
    inf=argv[0]
    stepi=int(argv[1])
    stepf=int(argv[2])
    steps=int(argv[3])
    ldict,datadict=read_log_lines(inf)
    neb_inps=get_neb_ordered(datadict,stepi,stepf,steps)
    for i,li in enumerate(neb_inps):
        with open (f'neb_in.{str(i)}','w') as f:
            for l in li:
                f.write(l)
    write_initial_coords(datadict[str(stepi)],ldict[str(stepi)])
    
if __name__ == "__main__":
    inp=sys.argv
    if len(inp) == 5:
        main(inp[1:])
    else:
        print("Improper inputs. Syntax should be:\n"
              "md2neb.py inputfile stepi stepf (stepsize)")
