from pymatgen.io.vasp.sets import MPRelaxSet
from fireworks import Firework
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet, ModifyIncar, WriteVaspStaticFromPrev
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from prlworkflows.input_sets import PRLRelaxSet, PRLStaticSet
from prlworkflows.prl_firetasks import WriteVaspFromIOSetPrevStructure, UpdateDisplacementDictForces, GeneratePhononDetour

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
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set))
        else:
            vasp_input_set = vasp_input_set or PRLRelaxSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type, gzip_output=False))
        t.append(PassCalcLocs(name=name))
        if db_insert:
            t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name, "metadata": metadata}))
        super(PRLOptimizeFW, self).__init__(t, parents=parents, name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)



class PRLStaticFW(Firework):
    def __init__(self, structure, name="static", vasp_input_set=None, dos=True, vasp_cmd="vasp", metadata=None,
                 prev_calc_loc=True, db_file=None, parents=None, phonon_detour=False, phonon_kwargs=None, **kwargs):
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
            Input set to use (for jobs w/no parents). Defaults to MPStaticSet() if None.
        dos : bool
            Whether to parse the electronic density of states. Defaults to True.
        vasp_cmd : str
            Command to run vasp.
        prev_calc_loc : (bool or str)
            If true (default), copies outputs from previous calc. If
            a str value, grabs a previous calculation output by name. If False/None, will create
            new static calculation using the provided structure.
        db_file (str): Path to file specifying db credentials.
        parents (Firework): Parents of this particular Firework. FW or list of FWS.
        prepare_phonon : bool
            If True, will add a task to create a phonon workflow as a detour.
        phonon_kwargs : dict
            Dictionary of arguments to pass to GeneratePhononDetour
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
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set))
        else:
            vasp_input_set = vasp_input_set or PRLStaticSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, parse_dos=True, additional_fields={"task_label": name, "metadata": metadata},))
        if phonon_detour:
            t.append(GeneratePhononDetour(**phonon_kwargs))
        super(PRLStaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class PRLPhononDisplacementFW(Firework):
    def __init__(self, structure, displacement_dict, name="phonon_displacement", vasp_input_set=None, smearing_type='gaussian', vasp_cmd="vasp",
                 parents=None, **kwargs):
        """
        Firework to calculate force sets.

        These always start from a fresh structure that has been displaced in a supercell.
        The forces calculated are pushed to the next FW using the spec.

        Parameters
        ----------
        structure : pymatgen.Structure
            Input structure.
        displacement_dict : dict
            A dictionary of
            ``{'direction': [1, 0, 0], 'displacement': np.array([0.01, 0, 0]), 'number': 0}``.
            This will be entered into the Firework's spec, updated with the forces
        name : str
            Name for the Firework.
        vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
            Input set to use. Defaults to PRLStaticSet() (with modifications) if None.
            The modifications are related to the ``smearing_type`` parameter. The input set
            must give accurate forces, (e.g. no tetrahedron smearing for metals).
        smearing_type : str
            This is only used when no ``vasp_input_set`` is passed. It must be one of 'gaussian',
            'methfessel-paxton', or 'tetrahedron'. The default is 'gaussian', which uses
            a SIGMA of 0.05. Using 'tetrahedron' gives a SIGMA of 0.05 and 'methfessel-paxton'
            a SIGMA of 0.2. Any further customizations should use a custom input set.
        vasp_cmd : str
            Command to run vasp.
        parents (Firework): Parents of this particular Firework. FW or list of FWS.
        \*\*kwargs: Other kwargs that are passed to Firework.__init__.

        Notes
        -----
        It is critical that an input set that can correctly give forces is used. Gaussian
        work for everything in principle, but is sensitive to the SIGMA parameter used.
        The tetrahedron method works well for insulators or semiconductors (which have full bands
        and no partial occupancies), but not metals, and requires no tuning of SIGMA.
        Methfessel-Paxton smearing gives accurate forces in metals, but should not be used for
        insulators and semiconductors. See the VASP manual for more details.
        """
        # put the displacement dict into the spec
        spec = kwargs.pop('spec', {})
        spec['displacement_dict'] = displacement_dict

        if vasp_input_set is None:
            if smearing_type == 'gaussian':
                vasp_input_set = PRLStaticSet(structure, user_incar_settings={'ISMEAR': 0, 'SIGMA': 0.05})
            elif smearing_type == 'tetrahedron':
                vasp_input_set = PRLStaticSet(structure)
            elif smearing_type == 'methfessel-paxton':
                vasp_input_set = PRLStaticSet(structure, user_incar_settings={'ISMEAR': 1, 'SIGMA': 0.2})
            else:
                raise ValueError('No vasp_input_set was passed and the smearing_type "{}" is not one of "gaussian", "tetrahedron" or "methfessel-paxton".'.format(smearing_type))

        tasks = []
        tasks.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        tasks.append(UpdateDisplacementDictForces())

        super(PRLPhononDisplacementFW, self).__init__(tasks, parents=parents, spec=spec,
                                                      name="{}-{}".format(structure.composition.reduced_formula, name),
                                                      **kwargs)
