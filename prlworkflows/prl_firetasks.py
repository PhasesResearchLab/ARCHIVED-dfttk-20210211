"""
Custom Firetasks for prlworkflows
"""
from pymatgen import Structure
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import load_class
from prlworkflows.analysis.phonon import get_all_force_sets, get_f_vib_phonopy


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
        unitcell = Structure.from_dict(fw_spec['unitcell'])
        supercell_matrix = fw_spec['supercell_matrix']

        temperatures, f_vib, s_vib, cv_vib = get_f_vib_phonopy(unitcell, supercell_matrix, disp_dicts)
        thermal_props_dict = {
            'volume': unitcell.volume,
            'F_vib': f_vib,
            'CV_vib': cv_vib,
            'S_vib': s_vib,
            'temperatures': temperatures,

        }

        return FWAction(mod_spec=[{'_push': {'f_vib': thermal_props_dict}}])

