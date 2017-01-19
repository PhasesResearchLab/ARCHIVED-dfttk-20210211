"""Describes running a NEB calculation in VASP

The logic here should be exclusively for NEB calculations that are abstracted from the acutal
inflection detection implementation. Until targets for release are found (e.g. pymatgen or atomate),
the style sholuld all follow atomate.
"""

from FireWorks import Firework, FiretaskBase
from pymatgen.analysis.transition_state import NEBAnalysis

class NEBTODbTask(FiretaskBase):
    """Uses NEBAnalysis from pymatgen to parse the outputs from an NEB calculation to """
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        # create an instance of NEBAnalysis. from_dir, current directory should be fine, but can also take a calc_loc from optional params. Each inner dir 00, 01... should have
        #     an out car and POSCAR
        # get the as_dict representation of the data
        # see if setting up the spline would be appropriate for finding the inflection or entering in the DB.
        # maybe this task would not even be needed because
        pass



class NEBFW(Firework):
    """For running an NEB calculation starting from a Henkelman and Jónsson setup."""
    def __init__(self, start_structure, end_structure, n_images=5, name="nudged elastic band", parents=None, vasp_cmd="vasp", db_file=None, **kwargs):
        """
        Args:
            start_structure: (Structure) starting structure
            end_structure: (Structure) ending structure
            n_images: (Int) number of images to create

        Returns:
            A Firework
        """
        # if parents, get the first and second images and replace the start and end structures. How to handle passing in the structures if we know there are parents? Just pass None?
        # get the structures from the previous calculations # should this be at the workflow level or FW?
        # calculate the images
        structures = start_structure.interpolate(end_structure, nimages=n_images+1, interpolate_lattices=True)
        # write vasp input for the start image in the main directory
        # create directories 00, 01, .. for the n_images and write POSCARs
        # run VASP - direct or with custodian (optional, of course). Should be able to override an OptmizeFW with the incar parameters
        # parse output task using the pymatgen analysis
        # make a tasks list, t
        super(NEBFW, self).__init__(t, parents=parents, name="{}-{}".format(
            start_structure.composition.reduced_formula, name), **kwargs)
