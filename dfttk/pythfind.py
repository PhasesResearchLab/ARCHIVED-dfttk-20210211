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
import shutil
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
        self.qhamode = args.qhamode
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        self.vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        #items = vasp_db.db['phonon'].find(properties=['metadata.tag','F_vib', 'CV_vib', 'S_vib'])
        #items = vasp_db.db['phonon'].find({F_vib: {$gt: 0}})
        #items = vasp_db.db['phonon'].find({'metadata.tag': "djdjd"})
        #items = vasp_db.db['phonon'].find({})
        self.items = (self.vasp_db).db[args.qhamode].find({})
        self.tags = []
        self.within = []
        self.containall = []
        self.containany = []
        if args.within is not None: self.within, tmp = formula2composition(args.within)
        if args.containall is not None: self.containall, tmp = formula2composition(args.containall)
        if args.containany is not None: self.containany, tmp = formula2composition(args.containany)

    def skipby(self, phase):
        els,tmp = formula2composition(phase.split('_')[0])
        if len (self.within) != 0:
            for e in els:
                if e not in self.within: return True
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
        for i in self.items:
            try:
                ii = len(i['S_vib'])
                mm = i['metadata']
                if ii > 0: 
                    if mm in hit:
                        count[hit.index(mm)] += 1
                        volumes[hit.index(mm)].append(i['volume'])
                    else:
                        hit.append(mm)
                        count.append(1)
                        volumes.append([i['volume']])
                        structure = Structure.from_dict(i['unitcell'])
                        formula_pretty = structure.composition.reduced_formula
                        sa = SpacegroupAnalyzer(structure)
                        phasename = formula_pretty+'_'\
                            + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                        phases.append(phasename)
            except:
               pass

        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        for i,m in enumerate(hit):
            static_calculations = (self.vasp_db).collection.find({'$and':[ {'metadata': m}, {'adopted': True} ]})
            skip_calc = False
            for calc in static_calculations:
                vol = calc['output']['structure']['lattice']['volume']
                if vol not in volumes[i] : 
                    skip_calc = True
                    break
            if skip_calc: continue
            if count[i] <= 5: continue
            if self.skipby(phases[i]): continue
            print (m, ":", phases[i], ": nV=", count[i])
            self.tags.append(m['tag'])


    def debye_find(self):
        hit = []
        phases = []
        count = []
        for i in self.items:
            try:
                ii = len(i['debye'])
                mm = i['metadata']
                if ii > 6: 
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
                        phases.append(phasename)
            except:
               pass
        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        for i,m in enumerate(hit):
            if self.skipby(phases[i]): continue
            print (m, ":", phases[i], count[i])
            self.tags.append(m['tag'])
