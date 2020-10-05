#!python

import os
from atomate.vasp.database import VaspCalcDb
from pymatgen.analysis.eos import Vinet, EOS
from pymatgen import Structure
from phonopy import Phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.phonopy import get_phonopy_structure
from dfttk.utils import sort_x_by_y
from dfttk.analysis.quasiharmonic import Quasiharmonic
import numpy as np
import json
import matplotlib.pyplot as plt
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn
        
class RetrieveDataFromDB(object):
    """docstring for RetrieveDataFromDB"""
    def __init__(self, db_file, properties, metadata, **kwargs):
        #super(RetrieveDataFromDB, self).__init__(db_file, properties, metadata, phonon)
        self.phonon = kwargs.get('phonon', False)
        self.properties = properties
        self.metadata = metadata
        self.t_min = kwargs.get('t_min', 5.)
        self.t_max = kwargs.get('t_max', 2000.)
        self.t_step = kwargs.get('t_step', 5.)
        self.poisson = kwargs.get('poisson', 0.363615)
        self.bp2gru = kwargs.get('bp2gru', 1)
        if phonon:
            qha_collectioin = 'qha_phonon'
        else:
            qha_collectioin = 'qha'
        property_map = {'qha': qha_collectioin, 'dos': 'tasks', 'eos': qha_collectioin, 
            'structure': qha_collectioin, 'phonon': 'phonon', 'debye': qha_collectioin}
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        self.vaspdb = vasp_db
        prop_items = list(vasp_db.db[property_map[properties]].find({'metadata': metadata}).sort('_id', 1))
        self.prop_items = prop_items

    def get_vol_energy_dos_from_static(self):
        static_calculations = self.vaspdb.db['tasks'].find({'$and':[ {'metadata': metadata}, {'adopted': True} ]})
        energies, volumes, dos_objs = [], [], []
        structure = None  # single Structure for QHA calculation
        for calc in static_calculations:
            energies.append(calc['output']['energy'])
            volumes.append(calc['output']['structure']['lattice']['volume'])
            dos_objs.append(self.vaspdb.get_dos(calc['task_id']))
            # get a Structure. We only need one for the masses and number of atoms in the unit cell.
            if structure is None:
                structure = Structure.from_dict(calc['output']['structure'])
        energies = sort_x_by_y(energies, volumes)
        dos_objs = sort_x_by_y(dos_objs, volumes)
        #structure = sort_x_by_y(structure)
        volumes = sorted(volumes)
        return (volumes, energies, dos_objs, structure)
        
    def get_eos(self, eos):
        self.head = ['volume', 'energy', 'energy_error']
        self.unit = ['A^3', 'eV', 'eV']
        self.data = np.vstack((eos.pop('volumes'), eos.pop('energies'), eos['error'].pop('difference'))).T
        self.parameter = eos

    def get_eos_from_static(self):
        self.head = ['volume', 'energy', 'energy_error']
        self.unit = ['A^3', 'eV', 'eV']
        (volumes, energies, dos_objs, structure) = self.get_vol_energy_dos_from_static()
        
        self.structure = structure
        eos = Vinet(volumes, energies)
        eos.fit()
        errors = eos.func(volumes) - energies
        self.data = np.vstack((volumes, energies, errors)).T
        self.parameter = {"b0_GPa": float(eos.b0_GPa), "b0": float(eos.b0), "b1": float(eos.b1),
                          "eq_volume": float(eos.v0),  "eq_energy": float(eos.e0), "name": "Vinet",
                          "error": {"sum_square_error": float(np.sum(np.square(errors)))}}
        self.volume = float(eos.v0)


    def get_qha(self, qha):
        self.head = ['Temperature', 'Optimum_volume', 'Gibbs_energy']
        self.unit = ['K', 'A^3', 'eV']
        self.data = np.vstack((qha.pop('temperatures'), qha.pop('optimum_volumes'), qha.pop('gibbs_free_energy'))).T
        self.parameter = qha

    def get_qha_from_static(self):
        self.head = ['Temperature', 'Optimum_volume', 'Gibbs_energy']
        self.unit = ['K', 'A^3', 'eV']
        
        (volumes, energies, dos_objs, structure) = self.get_vol_energy_dos_from_static()
        self.structure = structure

        f_vib = None

        # phonon properties
        if self.phonon:
            # get the vibrational properties from the FW spec
            phonon_calculations = list(self.vaspdb.db['phonon'].find({'metadata': self.metadata}).sort('_id', 1))
            vol_vol = [calc['volume'] for calc in phonon_calculations]  # these are just used for sorting and will be thrown away
            vol_f_vib = [calc['F_vib'] for calc in phonon_calculations]
            # sort them order of the unit cell volumes
            vol_f_vib = sort_x_by_y(vol_f_vib, vol_vol)
            vol_vol = sorted(vol_vol)
            f_vib = np.vstack(vol_f_vib)
            if len(volumes) > len(vol_vol):
                for i, voli in enumerate(vol_vol):
                    if abs(volumes[i] - voli) > 1e-3:
                        break
                if i == len(vol_vol)-1:
                    i = -1
                volumes.pop(index=i)
                energies.pop(index=i)

        qha_result = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=f_vib,
                            t_min=self.t_min, t_max=self.t_max, t_step=self.t_step,
                            poisson=self.poisson, bp2gru=self.bp2gru)

        eos = Vinet(volumes, energies)
        eos.fit()
        self.volume = float(eos.v0)

        qha = qha_result.get_summary_dict()
        qha['temperatures'] = qha['temperatures'].tolist()
        if not self.phonon:
            qha['poisson'] = self.poisson
            qha['bp2gru'] = self.bp2gru        

        self.data = np.vstack((qha.pop('temperatures'), qha.pop('optimum_volumes'), qha.pop('gibbs_free_energy'))).T
        self.parameter = qha        

    def write_output(self, **kwargs):
        if not self.prop_items:
            if self.properties == 'eos':
                self.get_eos_from_static()
            elif self.properties == 'qha':
                self.get_qha_from_static()
            else:
                raise ValueError
            self.formula = self.structure.composition.reduced_formula
            data_fn = "{}-{}-Vol{:.2f}.txt".format(self.formula, self.properties, self.volume)
            param_fn = "{}-{}-Vol{:.2f}.json".format(self.formula, self.properties, self.volume)
            with open(param_fn, "w+") as fjson:
                json.dump(self.parameter, fjson, indent=4)
            #write data
            with open(data_fn, "w+") as fdata:
                fdata.write("\t".join(self.head) + "\n")
                fdata.write("\t".join(self.unit) + "\n")
                for datai in self.data:
                    fdata.write("\t".join([str(item) for item in datai]) + "\n")
        else:
            for item in self.prop_items:
                if self.properties == 'phonon':
                    self.structure = Structure.from_dict(item['unitcell'])
                elif self.properties == 'dos':
                    self.structure = Structure.from_dict(item['output']['structure'])
                else:
                    self.structure = Structure.from_dict(item['structure'])
                self.formula = self.structure.composition.reduced_formula
                self.volume = self.structure.volume
                print(self.volume)
                data_fn = "{}-{}-Vol{:.2f}.txt".format(self.formula, self.properties, self.volume)
                param_fn = "{}-{}-Vol{:.2f}.json".format(self.formula, self.properties, self.volume)
                
                if self.properties == 'eos':
                    self.get_eos(item[self.properties])
                elif self.properties == 'qha':
                    if self.phonon:
                        self.get_qha(item['phonon'])
                    else:
                        self.get_qha(item['debye'])
                elif self.properties == 'phonon':
                    self.get_phonon(item, **kwargs)
                elif self.properties == 'debye':
                    self.get_debye(item['debye'])
                elif self.properties == 'dos':
                    self.get_dos(item)
                #'''
                #write json
                with open(param_fn, "w+") as fjson:
                    json.dump(self.parameter, fjson, indent=4)
                #write data
                with open(data_fn, "w+") as fdata:
                    fdata.write("\t".join(self.head) + "\n")
                    fdata.write("\t".join(self.unit) + "\n")
                    for datai in self.data:
                        fdata.write("\t".join([str(item) for item in datai]) + "\n")
                #'''

    def get_phonon(self, phonon, **kwargs):
        flag_savefig = kwargs.get('savefig', False)
        flag_savedata = kwargs.get('savedata', False)
        flag_band = kwargs.get('band', False)
        flag_dos = kwargs.get('phonon_dos', False)
        flag_pdos = kwargs.get('phonon_pdos', False)

        mesh = kwargs.get('mesh', [50, 50, 50])
        labels = kwargs.get('labels', None)

        if 'unitcell' not in phonon:
            raise FileNotFoundError('There is no phonon result. Please run phonon first.')
        self.head = ['Temperature', 'F_vib', 'CV_vib', 'S_vib']
        self.unit = ['K', 'eV', 'eV/K', 'eV/K']
        self.data = np.vstack((phonon.pop('temperatures'), phonon.pop('F_vib'), phonon.pop('CV_vib'), phonon.pop('S_vib'))).T
        filename = '{}-phonon-Vol{:.2f}'.format(self.formula, self.volume)

        unitcell = get_phonopy_structure(self.structure)
        supercell_matrix = phonon['supercell_matrix']
        force_constants = phonon['force_constants']
        ph = Phonopy(unitcell, supercell_matrix)
        ph.set_force_constants(force_constants)

        #for band structure
        if flag_band:
            if 'path' in kwargs:
                qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
                ph.run_band_structure(qpoints, path_connections=connections, labels=labels)
            else:
                ph.auto_band_structure()
            if flag_savefig:
                fig_band = ph.plot_band_structure()
                fig_band.savefig(fname = '{}-band.png'.format(filename))
                fig_band.close()
            if flag_savedata:
                ph.write_yaml_band_structure(filename = '{}-band.yaml'.format(filename))

        #for dos
        if flag_dos:
            ph.run_mesh(mesh)
            ph.run_total_dos()
            #phonon_dos_tmp = np.vstack((ph._total_dos._frequency_points, ph._total_dos._dos))
            #print(phonon_dos_tmp)
            #print(type(phonon_dos_tmp))
            if flag_savefig:
                fig_dos = ph.plot_total_dos()
                fig_dos.savefig(fname = '{}-dos.png'.format(filename))
                fig_dos.close()
            if flag_savedata:
                ph.write_total_dos(filename = '{}-dos.dat'.format(filename))
        #for pdos.
        if flag_pdos:
            ph.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
            ph.run_projected_dos()
            if flag_savefig:
                ph.plot_projected_dos().savefig(fname = '{}-pdos.png'.format(filename))
            if flag_savedata:
                ph.write_projected_dos(filename = '{}-pdos.dat'.format(filename))

        phonon.pop('_id')
        self.parameter = phonon

    def get_dos(self, dos):
        task_id = dos['task_id']
        dos_obj = self.vaspdb.get_dos(task_id)
        energies = dos_obj.energies
        efermi = dos_obj.efermi
        gap = dos_obj.get_gap()
        densities_up = dos_obj.densities[Spin.up]
        densities_dn = dos_obj.densities[Spin.down]
        #densities = densities_up + densities_dn
        densities = dos_obj.get_densities()
        #print('Fermi energy: {}'.format(efermi))
        #print('gap: {}'.format(gap))
        tdos = np.vstack((energies, densities))
        #print(tdos)
        #print(type(densities))

        plotter = DosPlotter()

        # Adds a DOS with a label.
        #plotter.add_dos("Total DOS", dos)

        # Alternatively, you can add a dict of DOSs. This is the typical
        # form returned by CompleteDos.get_spd/element/others_dos().
        #plotter.add_dos_dict({"dos1": dos1, "dos2": dos2})
        dosplot = plotter.add_dos_dict(dos_obj.get_element_dos())
        #dosplot = plotter.add_dos_dict(dos_obj.get_spd_dos())
        #plotter.show()

        self.head = ['energies', 'total_dos', 'total_dos_spin_up', 'total_dos_spin_down']
        self.unit = ['eV', '', '', '']
        self.data = np.vstack((energies, densities, densities_up, -densities_dn)).T
        self.parameter = dos['input']['incar']

    def get_debye(self, debye):
        self.head = ['Temperature', 'Optimum_volume', 'Gibbs_energy']
        self.unit = ['K', 'A^3', 'eV']
        self.data = np.vstack((debye.pop('temperatures'), debye.pop('optimum_volumes'), debye.pop('gibbs_free_energy'))).T
        self.parameter = debye


db_file = 'db.json'
#db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
properties = 'phonon'
#for Ni
metadata = {'tag': '217fecb5-4b27-4a7d-8ba9-1acc0bf5ffa9'}

phonon = True
band = True
phonon_dos = True

savefig = True
savedata = True

t_min = 5.
t_max = 2000.
t_step = 5.
poisson = 0.363615
bp2gru = 1
SETTINGS = {'phonon': phonon, 't_min': t_min, 't_max': t_max, 't_step': t_step, 'poisson': poisson, 'bp2gru': bp2gru}
data = RetrieveDataFromDB(db_file, properties, metadata, **SETTINGS)
data.write_output(band=band, phonon_dos=phonon_dos, savefig=savefig, savedata=savedata)