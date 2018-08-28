"""
Custom Firetasks for the DFTTK
"""
from pymatgen import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.eos import EOS
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import load_class, env_chk
from atomate.vasp.database import VaspCalcDb
from dfttk.analysis.phonon import get_f_vib_phonopy
from dfttk.analysis.quasiharmonic import Quasiharmonic
from dfttk.utils import sort_x_by_y
import numpy as np


@explicit_serialize
class WriteVaspFromIOSetPrevStructure(FiretaskBase):
    """
    Create VASP input files using implementations of pymatgen's VaspInputSet, overriding the
    Structure using a POSCAR in the current directory. An input set
    can be provided as an object or as a String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet object or a string
            name for the VASP input set (e.g., "StaticSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set, use this as a dict
            to specify kwargs for instantiating the input set parameters. For example, if you want
            to change the user_incar_settings, you should provide: {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object representation of a VaspInputSet
            rather than a String.
    """

    required_params = ["vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        struct = Structure.from_file('POSCAR')
        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']
            vis.structure = struct
        # if VaspInputSet String + parameters was provided
        else:
            vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])
            vis = vis_cls(struct, **self.get("vasp_input_params", {}))
        vis.write_input(".")


@explicit_serialize
class SupercellTransformation(FiretaskBase):
    """
    Transform a unitcell POSCAR to a supercell. Make a copy of the unitcell as a POSCAR-unitcell.

    This requires that a POSCAR is present in the current directory.
    """

    required_params = ['supercell_matrix']
    def run_task(self, fw_spec):
        unitcell = Structure.from_file('POSCAR')

        # create the unitcell file backup
        unitcell.to(filename='POSCAR-unitcell')

        # make the supercell and write to file
        unitcell.make_supercell(self['supercell_matrix'])
        unitcell.to(filename='POSCAR')


@explicit_serialize
class CalculatePhononThermalProperties(FiretaskBase):
    """
    Phonon-related Firetask to calculate force constants and F_vib.

    This requires that a vasprun.xml from a force constants run and
    a POSCAR-unitcell be present in the current directory.
    """

    required_params = ['supercell_matrix', 't_min', 't_max', 't_step', 'db_file', 'tag']
    optional_params = ['metadata']

    def run_task(self, fw_spec):

        tag = self["tag"]
        metadata = self.get('metadata', {})
        metadata['tag'] = tag

        unitcell = Structure.from_file('POSCAR-unitcell')
        supercell_matrix = self['supercell_matrix']
        temperatures, f_vib, s_vib, cv_vib, force_constants = get_f_vib_phonopy(unitcell, supercell_matrix, vasprun_path='vasprun.xml', t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'])
        if isinstance(supercell_matrix, np.ndarray):
            supercell_matrix = supercell_matrix.tolist()  # make serializable
        thermal_props_dict = {
            'volume': unitcell.volume,
            'F_vib': f_vib.tolist(),
            'CV_vib': cv_vib.tolist(),
            'S_vib': s_vib.tolist(),
            'temperatures': temperatures.tolist(),
            'force_constants': force_constants.tolist(),
            'metadata': metadata,
            'unitcell': unitcell.as_dict(),
            'supercell_matrix': supercell_matrix,
        }

        # insert into database
        db_file = env_chk(self["db_file"], fw_spec)
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        vasp_db.db['phonon'].insert_one(thermal_props_dict)



@explicit_serialize
class QHAAnalysis(FiretaskBase):
    """
    Do the quasiharmonic calculation from either phonon or Debye.

    Required params
    ---------------
    tag : str
        Tag to search the database for static calculations (energies, volumes, eDOS) from this job.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    phonon : bool
        True if f_vib comes from phonon calculations (in the spec). If False, it is calculated by the Debye model.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)

    Optional params
    ---------------
    poisson : float
        Poisson ratio, defaults to 0.25. Only used in Debye
    bp2gru : float
        Debye model fitting parameter for dBdP in the Gruneisen parameter. 2/3 is the high temperature
        value and 1 is the low temperature value. Defaults to 1.

    Notes
    -----
    Heterogeneity in the sources of E_0/F_el and F_vib is solved by sorting them according to increasing volume.
    """

    required_params = ["phonon", "db_file", "t_min", "t_max", "t_step", "tag"]

    optional_params = ["poisson", "bp2gru", "metadata"]

    def run_task(self, fw_spec):
        # handle arguments and database setup
        db_file = env_chk(self.get("db_file"), fw_spec)
        tag = self["tag"]

        vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)

        # get the energies, volumes and DOS objects by searching for the tag
        static_calculations = vasp_db.collection.find({"metadata.tag": tag})

        energies = []
        volumes = []
        dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
        structure = None  # single Structure for QHA calculation
        for calc in static_calculations:
            energies.append(calc['output']['energy'])
            volumes.append(calc['output']['structure']['lattice']['volume'])
            dos_objs.append(vasp_db.get_dos(calc['task_id']))
            # get a Structure. We only need one for the masses and number of atoms in the unit cell.
            if structure is None:
                structure = Structure.from_dict(calc['output']['structure'])

        # sort everything in volume order
        # note that we are doing volume last because it is the thing we are sorting by!

        energies = sort_x_by_y(energies, volumes)
        dos_objs = sort_x_by_y(dos_objs, volumes)
        volumes = sorted(volumes)


        qha_result = {}
        qha_result['structure'] = structure.as_dict()
        qha_result['formula_pretty'] = structure.composition.reduced_formula
        qha_result['metadata'] = self.get('metadata', {})
        qha_result['has_phonon'] = self['phonon']

        # phonon properties
        if self['phonon']:
            # get the vibrational properties from the FW spec
            phonon_calculations = list(vasp_db.db['phonon'].find({'metadata.tag': tag}))
            vol_vol = [calc['volume'] for calc in phonon_calculations]  # these are just used for sorting and will be thrown away
            vol_f_vib = [calc['F_vib'] for calc in phonon_calculations]
            # sort them order of the unit cell volumes
            vol_f_vib = sort_x_by_y(vol_f_vib, vol_vol)
            f_vib = np.vstack(vol_f_vib)
            qha = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=f_vib,
                                t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'],
                                poisson=self.get('poisson', 0.25), bp2gru=self.get('bp2gru', 1))
            qha_result['phonon'] = qha.get_summary_dict()
            qha_result['phonon']['temperatures'] = qha_result['phonon']['temperatures'].tolist()

        # calculate the Debye model results no matter what
        qha_debye = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=None,
                                  t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'],
                                  poisson=self.get('poisson', 0.25), bp2gru=self.get('bp2gru', 1))


        qha_result['debye'] = qha_debye.get_summary_dict()
        qha_result['debye']['temperatures'] = qha_result['debye']['temperatures'].tolist()

        # write to JSON for debugging purposes
        import json
        with open('qha_summary.json', 'w') as fp:
            json.dump(qha_result, fp)

        vasp_db.db['qha'].insert_one(qha_result)


@explicit_serialize
class EOSAnalysis(FiretaskBase):
    """
    Fit results from an E-V curve and enter the results the database

    Required params
    ---------------
    eos : str
        String name of the equation of state to use. See ``pymatgen.analysis.eos.EOS.MODELS`` for options.
    tag : str
        Tag to search the database for the volumetric calculations (energies, volumes) from this job.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.

    Optional params
    ---------------
    metadata : dict
        Metadata about this workflow. Defaults to an empty dictionary

    """

    required_params = ["eos", "db_file", "tag"]
    optional_params = ["metadata", ]

    def run_task(self, fw_spec):
        db_file = env_chk(self.get("db_file"), fw_spec)
        tag = self["tag"]
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        static_calculations = vasp_db.collection.find({"metadata.tag": tag})

        energies = []
        volumes = []
        structure = None  # single Structure for QHA calculation
        for calc in static_calculations:
            energies.append(calc['output']['energy'])
            volumes.append(calc['output']['structure']['lattice']['volume'])
            if structure is None:
                structure = Structure.from_dict(calc['output']['structure'])

        eos = EOS(self.get('eos'))
        ev_eos_fit = eos.fit(volumes, energies)
        equil_volume = ev_eos_fit.v0

        structure.scale_lattice(equil_volume)

        analysis_result = ev_eos_fit.results
        analysis_result['b0_GPa'] = float(ev_eos_fit.b0_GPa)
        analysis_result['structure'] = structure.as_dict()
        analysis_result['formula_pretty'] = structure.composition.reduced_formula
        analysis_result['metadata'] = self.get('metadata', {})
        analysis_result['energies'] = energies
        analysis_result['volumes'] = volumes


        # write to JSON for debugging purposes
        import json
        with open('eos_summary.json', 'w') as fp:
            json.dump(analysis_result, fp)

        vasp_db.db['eos'].insert_one(analysis_result)
