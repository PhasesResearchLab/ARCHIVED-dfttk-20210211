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
import dfttk.pyphon as pyphon
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
        self.qhamode = args.qhamode
        if args.qhamode == 'debye' : self.qhamode = 'qha'
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        self.vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        #items = vasp_db.db['phonon'].find(properties=['metadata.tag','F_vib', 'CV_vib', 'S_vib'])
        #items = vasp_db.db['phonon'].find({F_vib: {$gt: 0}})
        #items = vasp_db.db['phonon'].find({'metadata.tag': "djdjd"})
        #items = vasp_db.db['phonon'].find({})
        self.items = (self.vasp_db).db[args.qhamode].find({})
        self.tags = []
        self._Yphon = []
        self.within = []
        self.containall = []
        self.containany = []
        self.nV = args.nV
        self.get = args.get
        self.supercellN = args.supercellN
        self.toyphon = args.toyphon
        self.pyphon = args.pyphon
        self.t0 = args.t0
        self.t1 = args.t1
        self.td = args.td
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
        return self.tags, self._Yphon

    def toYphon(self,item, phase):
        volumes = []
        F_vib = []
        CV_vib = []
        S_vib = []
        T_vib = pythelec.T_remesh(self.t0, self.t1, self.td, 129)
        print ("extract the superfij.out used by Yphon ...")
        if not os.path.exists(phase):
            os.mkdir(phase)
        phdir = phase+'/Yphon'
        if not os.path.exists(phdir):
            os.mkdir(phdir)
        for i in (self.vasp_db).db[self.qhamode].find({'metadata.tag': item['metadata']['tag']}):
            if float(i['volume']) in volumes: continue
            volumes.append(float(i['volume']))
            structure = Structure.from_dict(i['unitcell'])
            natom = len(structure.sites)
            poscar = structure.to(fmt="poscar")
            unitcell_l = str(poscar).split('\n')
            supercell_matrix = i['supercell_matrix']
            supercell_structure = copy.deepcopy(structure)
            supercell_structure.make_supercell(supercell_matrix)
            natoms = len(supercell_structure.sites)
            poscar = supercell_structure.to(fmt="poscar")
            supercell_l = str(poscar).split('\n')
            vol = 'V{:010.6f}'.format(float(i['volume']))
            voldir = phdir+'/'+vol
            if not os.path.exists(voldir):
               os.mkdir(voldir)
            structure.to(filename=voldir+'/POSCAR')
            with open (voldir+'/metadata.json','w') as out:
                mm = item['metadata']
                mm['volume'] = i['volume']
                
                out.write('{}\n'.format(mm))
            with open (voldir+'/superfij.out','w') as out:
                for line in range (2,5):
                    out.write('{}\n'.format(unitcell_l[line]))
                for line in range (2,5):
                    out.write('{}\n'.format(supercell_l[line]))
                out.write('{} {}\n'.format(natoms, natoms//natom))
                for line in range (7,natoms+8):
                    out.write('{}\n'.format(supercell_l[line]))
                force_constant_matrix = np.array(i['force_constants'])
                hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)
                for ii in range(natoms):
                    for jj in range(natoms):
                        for x in range(3):
                            for y in range(3):
                                hessian_matrix[ii*3+x, jj*3+y] = -force_constant_matrix[ii,jj,x,y]
                for xx in range(natoms*3):
                    for yy in range(natoms*3-1):
                        out.write('{} '.format(hessian_matrix[xx,yy]))
                    out.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))
            if self.pyphon and self.get:
                cwd = os.getcwd()
                os.chdir( voldir )
                if not os.path.exists('vdos.out'):
                    cmd = "Yphon -tranI 2 -DebCut 0.5 " + " <superfij.out"
                    print(cmd, " at ", voldir)
                    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        universal_newlines=True)
                try:
                    if len(F_vib)==0:
                        print ("Calling yphon to get f_vib, s_vib, cv_vib at ", phdir)
                    with open("vdos.out", "r") as fp:
                        f_vib, U_ph, s_vib, cv_vib, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT \
                            = pyphon.vibrational_contributions(T_vib, dos_input=fp, energyunit='eV')
                except:
                    print ("Calling pyphon failed at ", voldir+'/vdos.out')
                    sys.exit()
                F_vib.append(f_vib)
                S_vib.append(s_vib)
                CV_vib.append(cv_vib)
                os.chdir( cwd )
        return {'T_vib':list(T_vib), 'volumes':volumes, 'F_vib':F_vib, 'S_vib':S_vib, 'CV_vib':CV_vib}

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
            for calc in static_calculations:
                vol = calc['output']['structure']['lattice']['volume']
                nS += 1
            if count[i] < self.nV: continue
            if self.supercellsize[i] < self.supercellN: continue
            sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, supercellsize: {:>3}, {}\n'.format(m, count[i], nS, self.supercellsize[i], phases[i]))
            if count[i]<6: continue
            if self.toyphon: 
                self.tags.append(m['tag'])
                self._Yphon.append(self.toYphon(ITEMS[i],phases[i]))
            else:
                self.tags.append(m['tag'])


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
            self.tags.append(m['tag'])
