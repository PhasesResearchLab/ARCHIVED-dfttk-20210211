#!/usr/bin/env python
from atomate.vasp.workflows.presets.core import wf_gibbs_free_energy
from atomate.vasp.powerups import add_modify_incar, remove_custodian
import numpy as np
from prlworkflows.utils import mp_structures_from_ids, get_launchpad

################################################################################
#                                CONFIGURATION                                 #
################################################################################

# Required configuration
use_api_interface = True  # set to True if you want to use the Materials Project API, False will load the POSCAR in the local dir
mp_structure_id = 'mp-134'  # enter the structure ID. Example: Al is 'mp-134'. See prlworkflows.utils for other options to get Structures.
use_phonopy = True
use_debye = False
is_conductor = True # TODO: needed with custodian?

# Optional configuration
launchpad_file_path = None  # absolute path to LaunchPad. If None, will load from FW_CONFIG_FILE variable
API_KEY = None  # Enter your API key here. If None, your PMG environment variable will be used

# calculation settings
custom_supercell_size = None  # if you want to use a custom supercell, pass a 3 element array, e.g. [3, 3, 3], otherwise use None. Passing None will calculate the minimum size required to get an 8x8x8 angstrom supercell.
deformations = None  # use custom deformations. Should be list of 3x3 transformation matrices (lists)
# for the INCAR settings, use a dictionary of {"setting": value}, e.g. {"EDIFFG": -0.05, "GGA": "PS"}
custom_incar_settings_all = None
custom_incar_settings_opt = None
custom_incar_settings_gibbs = None

# analysis settings
eos = "vinet"  # either "birch_murnaghan", "vinet" ...
t_max = 1000  # maximum temperature


################################################################################
#                                     RUN                                      #
################################################################################

def main():
    launchpad = get_launchpad(launchpad_file=launchpad_file_path)

    if use_api_interface:
        struct = mp_structures_from_ids([mp_structure_id], API_KEY=API_KEY)[0]
    else:
        from pymatgen import Structure

        struct = Structure.from_file('POSCAR')

    if custom_supercell_size:
        ss = custom_supercell_size
    else:
        lattice_sizes = [struct.lattice.a, struct.lattice.b, struct.lattice.c]
        ss = [np.ceil(8 / x) for x in lattice_sizes]
    struct.make_supercell(ss)

    c = {"t_max": t_max, "eos": eos}
    if deformations:
        c["deformations"] = deformations

    qhas = []
    if use_phonopy: qhas.append("phonopy")
    if use_debye: qhas.append("debye_model")
    if len(qhas) < 1:
        raise ValueError("Either use_phonopy or use_debye must be True.")

    for qha in qhas:
        c["qha_type"] = qha
        workflow = wf_gibbs_free_energy(struct, c)
        if is_conductor:
            workflow = add_modify_incar(workflow, modify_incar_params={"incar_update": {"SIGMA": 0.2, "ISMEAR": 1}})
        if custom_incar_settings_opt:
            workflow = add_modify_incar(workflow, modify_incar_params={"incar_update": custom_incar_settings_opt}, fw_name_constraint='optimization')
        if custom_incar_settings_all:
            workflow = add_modify_incar(workflow, modify_incar_params={"incar_update": custom_incar_settings_all})
        if custom_incar_settings_gibbs:
            workflow = add_modify_incar(workflow, modify_incar_params={"incar_update": custom_incar_settings_gibbs}, fw_name_constraint='gibbs')
        launchpad.add_wf(workflow)

if __name__ == '__main__':
    main()
