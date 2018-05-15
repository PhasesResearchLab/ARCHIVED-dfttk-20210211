"""
Custom Firetasks for prlworkflows
"""
from pymatgen import Structure
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import load_class, env_chk
from atomate.vasp.database import VaspCalcDb
from prlworkflows.analysis.phonon import get_all_force_sets, get_f_vib_phonopy
from prlworkflows.analysis.quasiharmonic import Quasiharmonic
from prlworkflows.utils import sort_x_by_y
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
            name for the VASP input set (e.g., "PRLStaticSet").

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
class UpdateDisplacementDictForces(FiretaskBase):
    """
    Phonon-related Firetask to get the forces from the vasprun.xml file in the current directory
    and pass this structure's displacement dict on to the next Firework.

    This requires that a displacement dictionary is present in the current Firework's spec.
    """
    def run_task(self, fw_spec):
        # get the forces, assuming a vasprun.xml is present in the current directory
        forces = get_all_force_sets(['vasprun.xml'])[0]

        disp_dict = fw_spec['displacement_dict']
        disp_dict['forces'] = forces

        return FWAction(mod_spec=[{'_push': {'displacement_dicts': disp_dict}}])


@explicit_serialize
class CalculatePhononThermalProperties(FiretaskBase):
    """
    Phonon-related Firetask to take displacement dictionaries, unitcell, and supercell matrix
    from the fw_spec and calculate F_vib.

    This requires that a list of ``displacement_dicts`` be present in the current Firework's spec.
    """
    def run_task(self, fw_spec):
        disp_dicts = fw_spec['displacement_dicts']
        # FireWorks unwraps forces arrays from NumPy to lists. We have to convert back otherwise we get errors in phonopy.
        force_sets = [np.array(ds.pop('forces')) for ds in disp_dicts]
        unitcell = fw_spec['unitcell']
        supercell_matrix = fw_spec['supercell_matrix']
        temperatures, f_vib, s_vib, cv_vib = get_f_vib_phonopy(unitcell, supercell_matrix, disp_dicts, force_sets=force_sets)
        thermal_props_dict = {
            'volume': unitcell.volume,
            'F_vib': f_vib,
            'CV_vib': cv_vib,
            'S_vib': s_vib,
            'temperatures': temperatures,
        }

        return FWAction(stored_data=thermal_props_dict, mod_spec=[{'_push': {'f_vib': thermal_props_dict}}])


@explicit_serialize
class GeneratePhononDetour(FiretaskBase):
    """
    A Firetask to create a phonon workflow from a single Firetask.

    Assumes the structure to start with is a POSCAR in the current directory.

    Required params
    ---------------
    supercell_matrix : numpy.ndarray
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]].

    Optional params
    ---------------
    smearing_type : str
        It must be one of 'gaussian', 'methfessel-paxton', or 'tetrahedron'. The default
        is 'methfessel-paxton', which uses a SIGMA of 0.2 and is well suited for metals,
        but it should not be used for semiconductors or insulators. Using 'tetrahedron'
        or 'gaussian' gives a SIGMA of 0.05. Any further customizations should use a custom workflow.
    displacement_distance : float
        Distance of each displacement. Defaults to 0.01, consistent with phonopy.
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.

    """

    required_params = ["supercell_matrix"]
    optional_params = ["smearing_type", "displacement_distance", "vasp_cmd", "t_step", "t_min", "t_max"]


    def run_task(self, fw_spec):
        from prlworkflows.prl_workflows import get_wf_phonon_single_volume

        smearing_type = self.get('smearing_type', 'methfessel-paxton')
        displacement_distance = self.get('displacement_distance', 0.01)
        vasp_cmd = self.get('vasp_cmd', None)

        phonon_wf = get_wf_phonon_single_volume(Structure.from_file('POSCAR'),
                                                self['supercell_matrix'],
                                                smearing_type=smearing_type,
                                                displacement_distance=displacement_distance,
                                                vasp_cmd=vasp_cmd,)

        return FWAction(detours=phonon_wf)


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

    optional_params = ["poisson", "bp2gru"]

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
        metadata = None  # single metadata for the tag to add to this job
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

        if self['phonon']:
            # get the vibrational properties from the FW spec
            vol_vol = [sp['volume'] for sp in fw_spec['f_vib']]  # these are just used for sorting and will be thrown away
            vol_f_vib = [sp['F_vib'] for sp in fw_spec['f_vib']]
            # sort them order of the unit cell volumes
            vol_f_vib = sort_x_by_y(vol_f_vib, vol_vol)
            f_vib = np.vstack(vol_f_vib)
        else:
            f_vib = None

        qha = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=f_vib,
                            t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'],
                            poisson=self.get('poisson', 0.25), bp2gru=self.get('bp2gru', 1))

        qha_summary = qha.get_summary_dict()
        qha_summary['temperatures'] = qha_summary['temperatures'].tolist()
        qha_summary['phonon'] = self['phonon']
        qha_summary['structure'] = structure.as_dict()
        qha_summary['formula_pretty'] = structure.composition.reduced_formula

        # write to JSON for debugging purposes
        import json
        with open('qha_summary.json', 'w') as fp:
            json.dump(qha_summary, fp)

        vasp_db.db['qha'].insert_one(qha_summary)


