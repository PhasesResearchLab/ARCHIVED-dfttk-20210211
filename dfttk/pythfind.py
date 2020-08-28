# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen import MPRester, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Potcar
#from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import dfttk.pythelec as pythelec
from dfttk.pythelec import thelecMDB
import warnings
import copy
import os
import sys
import subprocess
import shutil
import numpy as np
from fireworks.fw_config import config_to_dict
from atomate.vasp.database import VaspCalcDb
from dfttk.analysis.ywutils import formula2composition

class thfindMDB ():
    """
    find the metadata tag that has finished.

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
    """
    def __init__(self, args):
        if args.qhamode is not None:
            self.qhamode = args.qhamode
        else:
            self.qhamode = 'phonon'
        if args.qhamode == 'debye' : self.qhamode = 'qha'
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        self.vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        #items = vasp_db.db['phonon'].find(properties=['metadata.tag','F_vib', 'CV_vib', 'S_vib'])
        #items = vasp_db.db['phonon'].find({F_vib: {$gt: 0}})
        #items = vasp_db.db['phonon'].find({'metadata.tag': "djdjd"})
        #items = vasp_db.db['phonon'].find({})
        self.items = (self.vasp_db).db[self.qhamode].find({})
        self.tags = []
        self._Yphon = []
        self.within = []
        self.containall = []
        self.containany = []
        self.excludeall = []
        self.excludeany = []
        self.nV = args.nV
        self.get = args.get
        self.supercellN = args.supercellN
        self.t0 = args.t0
        self.t1 = args.t1
        self.td = args.td
        self.findbandgap = args.findbandgap
        if args.within is not None: self.within, tmp = formula2composition(args.within)
        if args.containall is not None: self.containall, tmp = formula2composition(args.containall)
        if args.containany is not None: self.containany, tmp = formula2composition(args.containany)
        if args.excludeall is not None: self.excludeall, tmp = formula2composition(args.excludeall)

    def skipby(self, phase):
        els,tmp = formula2composition(phase.split('_')[0])
        if len (self.within) != 0:
            for e in els:
                if e not in self.within: return True
        if len (self.excludeall) != 0:
            for e in self.excludeall:
                rx = e not in els
                if rx: break
            if not rx: return True
        if len (self.excludeany) != 0:
            for e in self.excludeany:
                return True
        if len (self.containall) != 0:
            for e in self.containall:
               if e not in els: return True
            return False
        if len (self.containany) != 0:
            for e in self.containany:
                if e in els: return False
            return True
        return False
    
    def run_console(self):
        if self.qhamode=='phonon': self.phonon_find()
        else: self.debye_find()
        return self.tags


    def phonon_find(self):
        hit = []
        count = []
        phases = []
        volumes = []
        ITEMS = []
        self.supercellsize = []
        for i in self.items:
            try:
                ii = len(i['S_vib'])
                mm = i['metadata']
            except:
                continue
            if ii <= 0: continue 
            if mm in hit:
                if i['volume'] not in volumes[hit.index(mm)]:
                    volumes[hit.index(mm)].append(i['volume'])
                    count[hit.index(mm)] += 1
            else:
                ITEMS.append(i)
                hit.append(mm)
                count.append(1)
                volumes.append([i['volume']])
                structure = Structure.from_dict(i['unitcell'])
                natoms = len(structure.sites)
                supercell_matrix = i['supercell_matrix']
                self.supercellsize.append(natoms*int(np.linalg.det(np.array(supercell_matrix))+.5))
                formula_pretty = structure.composition.reduced_formula
                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)

        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        for i,m in enumerate(hit):
            if self.skipby(phases[i]): continue
            static_calculations = (self.vasp_db).collection.\
                find({'$and':[ {'metadata.tag': m['tag']}, {'adopted': True} ]})
            nS = 0
            gapfound = False
            potsoc = ""
            for calc in static_calculations:
                vol = calc['output']['structure']['lattice']['volume']
                if potsoc=="":
                    pot = calc['input']['pseudo_potential']['functional'].upper()
                    try:
                        if calc['input']['incar']['LSORBIT']: potsoc = pot +"SOC"
                    except:
                        potsoc = pot
                    pname = phases[i].split('#')
                    if len(pname)>1: phases[i] = pname[0]+potsoc+'#'+pname[1]
                    else: phases[i] = pname[0]+potsoc
                nS += 1
                bandgap = calc['output']['bandgap']
                if not gapfound: gapfound = float(bandgap) > 0.0
            if self.findbandgap:
                #if gapfound: print ("eeeeee", gapfound, bandgap, phases[i])
                if gapfound: sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, supercellsize: {:>3}, {}\n'.format(m, count[i], nS, self.supercellsize[i], phases[i]))
            else:
                if count[i] < self.nV: continue
                if self.supercellsize[i] < self.supercellN: continue
                sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, supercellsize: {:>3}, {}\n'.format(m, count[i], nS, self.supercellsize[i], phases[i]))
                if count[i]>=6: self.tags.append({'tag':m['tag'],'phasename':phases[i]})


    def debye_find(self):
        hit = []
        phases = []
        count = []
        for i in self.items:
            try:
                ii = len(i['debye'])
                mm = i['metadata']
            except:
                continue
            if ii < 6: continue
            if mm in hit:
                count[hit.index(mm)] += 1
            else:
                hit.append(mm)
                count.append(1)
                structure = Structure.from_dict(i['structure'])
                formula_pretty = structure.composition.reduced_formula
                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)
        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        for i,m in enumerate(hit):
            if self.skipby(phases[i]): continue
            print (m, ":", phases[i])
            self.tags.append({'tag':m['tag'],'phasename':phases[i]})
