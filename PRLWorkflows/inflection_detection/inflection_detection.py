"""Implements an inflection detection workflow and customized tasks, per the van de Walle 2015 paper

The code here should only pertain to the specific implementation for the inflection detection method.
Any code that can be abstracted out or otherwise generalized to be useful should be refactored.
After this can become public, true abstraction and a preset workflow can be written for submission
to atomate.
"""

from FireWorks import Workflow, FiretaskBase, Firework
from PRLWorkflows.inflection_detection.neb import NEBFW

def calculate_inflection_energy(energies):
    """Use sympy to construct a fit where the inflection point can be calculated.

    Possibly NEBAnalysis spline could be used to construct the fit. Perhaps sympy can directly
    calculate the inflection from the list of points.

    Args:
        energies (list of float): energies to calculate the inflection

    Returns:
        A float of the inflection energy
    """
    pass

class InflectionDetectionToDbTask(FiretaskBase):
    """Finds the inflection point from energies that are in the database.

    Energies will come directly from database entries from required params. The resulting inflection
    point energy is entered into the database in a new collection, along with other relevant metadata.
    """
    required_params = []

    def run_task(self, fw_spec):
        # take the list of database entries from the required params. See how the GibbsFreeEnergyToDbTask works for this.
        # calculate the inflection energy
        # enter it in the DB with relevant metadata. Again see the Gibbs task.
        pass



def get_wf_inflection_detection(structure, n_images, vasp_cmd="vasp", db_file=None):
    """Get an inflection detection workflow

    Args:
        structure (Structure): unstable structure
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        A Workflow

    Example:
        >>> struct = Structure.from_file('POSCAR-FCC-W')
        >>> get_wf_inflection_detection(struct, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<") # use env_chk
        # returns the workflow
    """
    # fw0 = OptimizeFW() #customize with ISIF 7.
    # fw1 = OptimizeFW() # customize with ISIF 3
    # fw2 = NEBFW() # set parents = 0, 1. set n_images. Which structures to set for start and end? "None"?
    # fw3 to fw(3+n_images+2) = StaticFW for each image using a loop. Set the parents for each FW to 2
    # fw(3+n_images+3) = Firework(InflectionDetectionToDbTask(), name="Inflection Detection").
    # (above) Calculate the inflection based on the energies from static calculations. Parents are range(3..(3+n_images+3))
    pass