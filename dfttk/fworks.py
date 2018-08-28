import warnings
from uuid import uuid4

from fireworks import Firework
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from dfttk.input_sets import RelaxSet, StaticSet, ForceConstantsSet
from dfttk.ftasks import WriteVaspFromIOSetPrevStructure, SupercellTransformation, CalculatePhononThermalProperties


class OptimizeFW(Firework):
    """
    Optimize the given structure.

    Results are not entered into the database.

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
    def __init__(self, structure, name="structure optimization", vasp_input_set=None, job_type="normal",
                 vasp_cmd="vasp", metadata=None, override_default_vasp_params=None, db_file=None,
                 force_gamma=True, prev_calc_loc=True, parents=None, db_insert=False, **kwargs):

        metadata = metadata or {}
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or RelaxSet(structure, force_gamma=force_gamma,
                                                       **override_default_vasp_params)
        t = []
        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set))
        else:
            vasp_input_set = vasp_input_set or RelaxSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type, gzip_output=False))
        t.append(PassCalcLocs(name=name))
        if db_insert:
            t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name, "metadata": metadata}))
        super(OptimizeFW, self).__init__(t, parents=parents, name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)



class StaticFW(Firework):
    """
    Standard static calculation Firework - either from a previous location or from a structure.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to StaticSet() if None.
    vasp_cmd : str
        Command to run vasp.
    prev_calc_loc : (bool or str)
        If true (default), copies outputs from previous calc. If
        a str value, grabs a previous calculation output by name. If False/None, will create
        new static calculation using the provided structure.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    \*\*kwargs : dict
        Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, name="static", vasp_input_set=None, vasp_cmd="vasp", metadata=None,
                 prev_calc_loc=True, db_file=None, parents=None, **kwargs):

        # TODO: @computron - I really don't like how you need to set the structure even for
        # prev_calc_loc jobs. Sometimes it makes appending new FWs to an existing workflow
        # difficult. Maybe think about how to remove this need? -computron
        metadata = metadata or {}
        vasp_input_set = vasp_input_set or StaticSet(structure)

        t = []

        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set))
        else:
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, parse_dos=True, additional_fields={"task_label": name, "metadata": metadata},))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class PhononFW(Firework):
    """
    Calculation of phonon thermal properties by direct calculation of force constants.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    supercell_matrix:
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]].
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to ForceConstantsSet() if None.
    vasp_cmd : str
        Command to run vasp.
    prev_calc_loc : (bool or str)
        If true (default), copies outputs from previous calc. If
        a str value, grabs a previous calculation output by name. If False/None, will create
        new static calculation using the provided structure.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    \*\*kwargs : dict
        Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, supercell_matrix, t_min=5, t_max=2000, t_step=5,
                 name="phonon", vasp_input_set=None,
                 vasp_cmd="vasp", metadata=None, tag=None,
                 prev_calc_loc=True, db_file=None, parents=None,
                 **kwargs):

        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            warnings.warn('No ``tag`` was passed explicitly or in ``metadata`` to PhononFW. In order to find this Firework later, you should assign one. This was assigned: {}'.format(tag))
            metadata['tag'] = tag

        vasp_input_set = vasp_input_set or ForceConstantsSet(structure)

        t = []

        # We need to get the POSCAR from the previous run or from the passed Structure
        # so it can be transformed to a supercell in the next step
        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
        else:
            # write the input set first, just to get the POSCAR file in the directory
            # the other inputs will get overridden by WriteVaspFromIOSetPrevStructure
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(SupercellTransformation(supercell_matrix=supercell_matrix))
        t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(PassCalcLocs(name=name))
        t.append(CalculatePhononThermalProperties(supercell_matrix=supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag, metadata=metadata))

        super(PhononFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)
