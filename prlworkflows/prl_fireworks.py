from pymatgen.io.vasp.sets import MPRelaxSet
from fireworks import Firework
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet, ModifyIncar, WriteVaspStaticFromPrev
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from prlworkflows.input_sets import PRLRelaxSet
import warnings

class PRLOptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization", vasp_input_set=None, job_type="normal",
                 vasp_cmd="vasp", metadata=None, override_default_vasp_params=None, db_file=None,
                 force_gamma=True, prev_calc_loc=True, parents=None, db_insert=False, **kwargs):
        """
        Optimize the given structure.
        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, grabs a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            db_insert : bool
                Whether to insert the task into the database. Defaults to False.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        metadata = metadata or {}
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or PRLRelaxSet(structure, force_gamma=force_gamma,
                                                       **override_default_vasp_params)
        t = []
        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev())
        else:
            vasp_input_set = vasp_input_set or PRLRelaxSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type))
        t.append(PassCalcLocs(name=name))
        if db_insert:
            t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name, "metadata": metadata}))
        super(PRLOptimizeFW, self).__init__(t, parents=parents, name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)



class PRLStaticFW(Firework):
    def __init__(self, structure, name="static", vasp_input_set=None, vasp_cmd="vasp", metadata=None,
                 prev_calc_loc=True, db_file=None, parents=None, **kwargs):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure 
                is only used to set the name of the FW and any structure with the same composition 
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If 
                a str value, grabs a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        # TODO: @computron - I really don't like how you need to set the structure even for
        # prev_calc_loc jobs. Sometimes it makes appending new FWs to an existing workflow
        # difficult. Maybe think about how to remove this need? -computron
        metadata = metadata or {}

        t = []

        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev())
        else:
            vasp_input_set = vasp_input_set or MPStaticSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name, "metadata": metadata}))
        super(PRLStaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


