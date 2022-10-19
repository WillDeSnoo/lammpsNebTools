#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 09:00:14 2022

@author: willdesnoo
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def column(matrix, i):
    return [row[i] for row in matrix]


class NebOutData():
    def __init__(self):
        self.cNEBstart = None
        self.cNEBindex = None
        self.datadict = {}
        self.nreps = int

    def _list2dict(self, list2d):
        keys = list2d[0]
        for i, key in enumerate(keys):
            self.datadict[key] = list2d[1:, i].astype(float)

    def read_out_data(self, fname):
        dbool = False
        with open(fname) as out_f:
            out_lines = out_f.readlines()
            out_dat = []
            for i, l in enumerate(out_lines):
                if "Setting up climbing" in l:
                    self.cNEBstart = int(out_dat[-1][0])
                    self.cNEBindex = len(l)
                    dbool = False
                    continue
                if dbool:  # dbool=True if line contains data
                    spl = l.split()
                    out_dat.append(spl)
                if "Setting up regular NEB ..." in l:
                    dbool = True
                    continue

                elif not dbool:
                    if "Step MaxReplicaForce" in l:
                        dbool = True
                        continue

        # Fixing keys
        del out_dat[0][-3:]
        c = 3
        while len(out_dat[0]) < len(out_dat[1]):
            out_dat[0].append(f"RD{str(c)}")
            out_dat[0].append(f"PE{str(c)}")
            c += 1
        print(f"{str(c-1)} replicas detected.")
        self.nreps = c-1
        self._list2dict(np.asarray(out_dat))

    def plot_data(self, *iny, inx="", mintozero=False, ftol=None):
        ys = []
        if not inx:  # If x is not specified, "Step" is default x
            inx = "Step"
            
        absmin=min(self.datadict[str(iny[0])])
        for y in iny:
            currmin=min(self.datadict[str(y)])
            absmin = [currmin if currmin < absmin else absmin]
        for y in iny:
            if mintozero:
                ys.append([x-absmin
                          for x in self.datadict[str(y)]])
            else:
                ys.append(self.datadict[str(y)])
        x = self.datadict[inx]

        # Plotting setup
        # plt.ion()
        fig, ax = plt.subplots()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(int(max(x)/5)))
        plt.xticks(rotation='vertical')
        if ftol:
            plt.axhline(y=ftol, color='black', label='ftol', linestyle='--')
        if self.cNEBstart:
            plt.axvline(x=self.cNEBstart, color='black',
                        label='cNEB begins', linestyle='--')
        for i, y in enumerate(ys):
            ax.plot(x, y, label=f"PE{str(i+1)}")

        plt.xlabel(inx)
        if ("PE" in iny[0]):
            plt.ylabel("Potential Energy (kJ/mol)")
        elif ("Force" in iny[0]):
            plt.ylabel("Force (eV/Å)")
        plt.legend()

    def plot_finalE(self, pt='cneb'):
        PEs, RDs = [], []
        for i in range(1, self.nreps+1):
            PEs.append(self.datadict[f"PE{str(i)}"])
            RDs.append(self.datadict[f"RD{str(i)}"])
        fig, ax = plt.subplots()
        if pt == 'cneb': 
            plt.plot(column(RDs, -1), column(PEs, -1), marker='o')
            Ea=max(column(PEs,-1))-PEs[0][-1]
            print("cNEB Ea: {0:.2f} kJ/mol".format(Ea))
        elif pt == 'neb':
            plt.plot(column(RDs, self.cNEBindex), column(PEs, self.cNEBindex), marker='o')
            Ea=max(column(PEs,self.cNEBindex))-PEs[0][self.cNEBindex]
            print("NEB Ea: {0:.2f} kJ/mol".format(Ea))
            #plt.plot(column(RDs, 10), column(PEs, 10), marker='o')
        else:
            plt.plot(column(RDs, int(pt)), column(PEs, int(pt)), marker='o')
            
        plt.xlabel("Reaction Progress")
        plt.ylabel("Energy (kJ/mol")
        #plt.legend()
        

def main(argv):
    neb = NebOutData()
    neb.read_out_data(argv[0])
    dd = neb.datadict
    keys = dd.keys()
    nreps = neb.nreps
    print(keys)
    uin = input("Please input which key you would like to plot:\n"
                "Alternatively, type the following number for each listed option:\n"
                "0: Exit.\n"
                "1: Plot all potential energies.\n"
                "2: Plot Max Replica Force.\n"
                "3: Plot final replica energies\n"
                "4: Plot all options simultaneously")
    if uin == "0":
        print("Exiting...")
    elif uin == "1":
        PEs = []
        for i in range(1, nreps+1):
            PEs.append(f"PE{str(i)}")
        neb.plot_data(*PEs)
        plt.show()

    elif uin == "2":
        neb.plot_data("MaxReplicaForce", ftol=0.025)
        plt.show()

    elif uin == "3":
        neb.plot_finalE()
        plt.show()

    elif uin == "4":
        PEs = []
        for i in range(1, nreps+1):
            PEs.append(f"PE{str(i)}")
        neb.plot_data(*PEs, mintozero=False)
        neb.plot_data("MaxReplicaForce", ftol=0.025)
        if neb.cNEBstart:
            neb.plot_finalE(pt='neb')
        neb.plot_finalE()
        plt.show()

    elif uin not in keys:
        print("Not a valid key, please try again")

    else:
        print("Plotting...")
        neb.plot_data(uin)
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
