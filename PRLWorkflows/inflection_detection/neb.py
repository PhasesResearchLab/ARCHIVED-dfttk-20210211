"""Describes running a NEB calculation in VASP

The logic here should be exclusively for NEB calculations that are abstracted from the acutal
inflection detection implementation. Until targets for release are found (e.g. pymatgen or atomate),
the style sholuld all follow atomate.
"""

from FireWorks import Firework, FiretaskBase
from atomate.vasp.firetasks.write_inputs import WriteVaspFromPMGObjects, WriteVaspFromIOSet
from atomate.vasp.firetasks.run_calc import RunVaspDirect
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.analysis.transition_state import NEBAnalysis

class NEBTODbTask(FiretaskBase):
    """Uses NEBAnalysis from pymatgen to parse the outputs from an NEB calculation to a DB

    Optional params:
        base_calc_dir (str): path to dir (on current filesystem) that contains top level VASP
            output files. Default: use current working directory.
        calc_dirs ([str]): List of calculation directories for sub-calculations
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        defuse_unsuccessful (bool): Defuses children fireworks if VASP run state
            is not "successful"; i.e. both electronic and ionic convergence are reached.
            Defaults to True.
    """
    optional_params = ["db_file", "additional_fields", "base_calc_dir", "calc_dirs", "calc_loc", "defuse_unsuccesssful"]

    def run_task(self, fw_spec):
        # create an instance of NEBAnalysis. from_dir, current directory should be fine, but can also take a calc_loc from optional params. Each inner dir 00, 01... should have
        #     an out car and POSCAR
        # get the as_dict representation of the data
        # see if setting up the spline would be appropriate for finding the inflection or entering in the DB.
        # maybe this task would not even be needed because
        pass



class NEBFW(Firework):
    """For running an NEB calculation starting from a Henkelman and Jónsson setup."""
    def __init__(self, start_structure, end_structure, n_images=5, name="nudged elastic band", use_parent_structures=False, parents=None, vasp_cmd="vasp", db_file=None, **kwargs):
        """
        Args:
            start_structure (Structure): starting structure
            end_structure (Structure): ending structure

        Kwargs:
            n_images (Int): number of images to create
            use_parent_structures (Bool): if True, will use structures from 2 parents as start and end
                Default: False
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__
        Returns:
            A Firework
        """
        if use_parent_structures:
            if len(parents) != 2:
                raise ValueError("NEBFW must be 2 to use parent structures. {0} parents passed.".format(len(parents)))
            start_structure = parents[0].spec # TODO: get the output structure from a FW
            end_structure = parents[1].spec # TODO: get the output structure from a FW

        # calculate the images
        structures = start_structure.interpolate(end_structure, nimages=n_images+1, interpolate_lattices=True)

        t = []
        # write vasp input for the start image in the main directory
        # TODO: Do we need to turn on LNEBCELL or other flags? http://theory.cm.utexas.edu/vtsttools/neb.html#nudged-elastic-band-options
        uis = {"user_incar_settings":{"incar_update":{"ICHAIN": 0, "IMAGES":n_images}}}
        vasp_input_set = MPStaticSet(structures[0], override_default_vasp_params=uis)
        t.append(WriteVaspFromIOSet(structures[0]), vasp_input_set=vasp_input_set)
        # create directories 00, 01, .. for the n_images and write POSCARs
        calc_dirs = []
        for i, struct in enumerate(structures):
            path = '{:02d}'.format(i)
            calc_dirs.append(path)
            # TODO: how to add folder path to the write input task
            t.append(WriteVaspFromPMGObjects(poscar=Poscar(struct)))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(NEBTODbTask(db_file=db_file, additional_fields={"task_label": name, "n_images": n_images}, calc_dirs=calc_dirs))
        super(NEBFW, self).__init__(t, parents=parents, name="{}-{}".format(
            start_structure.composition.reduced_formula, name), **kwargs)
