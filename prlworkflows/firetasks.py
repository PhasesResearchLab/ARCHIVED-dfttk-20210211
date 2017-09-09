import os

from fireworks import explicit_serialize, FiretaskBase, FWAction
from pymatgen import Structure

from atomate.utils.utils import get_logger

logger = get_logger(__name__)

# TODO: Not safe for ssh. Use atomate.utils.fileio.FileClient
@explicit_serialize
class CheckSymmetry(FiretaskBase):
    """
    Checks the symmetry of the final structure against the starting structure.

    If the symmetry is the same, the ``pass_action`` will be invoked. The default is to continue the
    Firework and Workflow. If the symmetry has broken with respect to the tolerances, then the
    ``fail_action`` will be invoked, which defaults to exit and defuse the workflow by default.

    Required params:
        (none)
    Optional params:
        calc_dir (str): path containing the VASP input/output files. Defaults to current directory.
        symprec (float): Tolerance for symmetry finding. Defaults to 1e-2.
        angle_tolerance (float): Angle tolerance for symmetry finding. Defaults to 5 degrees.
        rename_input_on_fail (bool): Rename POSCAR -> CONTCAR on detection of broken symmetry.
            Useful for compatibility with ``CopyVaspOutputs``. Defaults to False.
        pass_action (FWAction): FWAction on success (symmetry preserved). Defaults to add stored data and continue.
        fail_action (FWAction): FWAction on failure (symmetry changed). Defaults to exit the FW and defuse the WF.
    """

    required_params = []
    optional_params = ["calc_dir", "symprec", "angle_tolerance", "rename_on_fail", "pass_action", "fail_action"]

    def run_task(self, fw_spec):
        calc_dir = self.get("calc_dir", ".")

        def find_file(filename):
            """Find file by filename with support for extensions."""
            possible_exts = ['', '.gz']
            for ext in possible_exts:
                if os.path.exists(filename+ext):
                    return filename+ext
            raise FileNotFoundError('File {} does not exist with any of the extensions: {}'.format(filename, possible_exts))

        # get initial structure, preferring POSCAR.orig to POSCAR
        try:
            input_file = find_file(os.path.join(calc_dir, 'POSCAR.orig'))
        except FileNotFoundError:
            input_file = (os.path.join(calc_dir, 'POSCAR'))
        initial_structure = Structure.from_file(input_file)
        contcar_file = find_file(os.path.join(calc_dir, 'CONTCAR'))
        final_structure = Structure.from_file(contcar_file)

        symprec = self.get("symprec", 1e-2)
        angle_tolerance = self.get("angle_tolerance", 5)
        # returns tuple of (sg symbol, number)
        initial_structure_sg_info = initial_structure.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
        final_structure_sg_info = final_structure.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)

        stored_data = {'initial_structure':
                           {'spacegroup_symmetry': initial_structure_sg_info[0],
                            'structure': initial_structure.as_dict()},
                       'final_structure':
                           {'spacegroup_symmetry': final_structure_sg_info[0],
                            'structure': final_structure.as_dict()}
                       }
        pass_action = self.get("pass_action", FWAction(stored_data=stored_data))
        fail_action = self.get("fail_action", FWAction(stored_data=stored_data, exit=True, defuse_workflow=True))

        if initial_structure_sg_info[0] != final_structure_sg_info[0]:
            logger.info("CheckSymmetry: symmetry of initial structure ({}) is different from the final structure ({}).".format(initial_structure_sg_info[0], final_structure_sg_info[0]))
            if self.get("rename_input_on_fail", False):
                os.rename(input_file, contcar_file)
                logger.info("CheckSymmetry: Moving {} to {}".format(input_file, contcar_file))
            return fail_action
        else:
            return pass_action
