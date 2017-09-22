from pymatgen.io.vasp.sets import MPRelaxSet
from fireworks import Firework
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.run_calc import RunVaspDirect

class OptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization", vasp_input_set=None,
                 vasp_cmd="vasp", isif=None, override_default_vasp_params=None, db_file=None,
                 force_gamma=True, parents=None, **kwargs):
        """
        Optimize the given structure.
        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            isif : int
                Shortcut to override the ISIF parameter. Defaults to None.
                Will take precedent over override_default_vasp_params
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        if isif:
            override_default_vasp_params.update({"ISIF": isif})
        vasp_input_set = vasp_input_set or MPRelaxSet(structure, force_gamma=force_gamma,
                                                      **override_default_vasp_params)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name),
                                         **kwargs)
