import six

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.sets import DictSet, get_vasprun_outcar, get_structure_from_prev_run, _load_yaml_config

class PRLRelaxSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['INCAR'].pop('ENCUT')
    CONFIG['KPOINTS'].update({
        'grid_density': 6000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].update({
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'PREC': 'HIGH',
        'ALGO': 'NORMAL'
    })

    def __init__(self, structure, **kwargs):
        self.kwargs = kwargs
        super(PRLRelaxSet, self).__init__(
            structure, PRLRelaxSet.CONFIG, **kwargs)

class PRLStaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'ENCUT': 520, # MP compatibility
        'ISMEAR': -5,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ENCUT': 520,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": True,
        "LORBIT": 11,
        "LVHAR": True,
        "LWAVE": False,
        "ICHARG": 0,
    })

    def __init__(self, structure, prev_incar=None, prev_kpoints=None,
                 lepsilon=False, lcalcpol=False, grid_density=8000, **kwargs):
        super(PRLStaticSet, self).__init__(structure, PRLStaticSet.CONFIG, **kwargs)

        if isinstance(prev_incar, six.string_types):
            prev_incar = Incar.from_file(prev_incar)
        if isinstance(prev_kpoints, six.string_types):
            prev_kpoints = Kpoints.from_file(prev_kpoints)

        self.prev_incar = prev_incar
        self.prev_kpoints = prev_kpoints
        self.grid_density = grid_density
        self.structure = structure
        self.kwargs = kwargs
        self.lepsilon = lepsilon
        self.lcalcpol = lcalcpol

    @property
    def incar(self):
        parent_incar = super(PRLStaticSet, self).incar
        incar = Incar(self.prev_incar) if self.prev_incar is not None else \
            Incar(parent_incar)

        incar.update(
            {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "ALGO": "Normal"})

        if self.lepsilon:
            incar["IBRION"] = 8
            incar["LEPSILON"] = True

            # LPEAD=T: numerical evaluation of overlap integral prevents
            # LRF_COMMUTATOR errors and can lead to better expt. agreement
            # but produces slightly different results
            incar["LPEAD"] = True

            # Note that DFPT calculations MUST unset NSW. NSW = 0 will fail
            # to output ionic.
            incar.pop("NSW", None)
            incar.pop("NPAR", None)

        if self.lcalcpol:
            incar["LCALCPOL"] = True

        for k in ["MAGMOM", "NUPDOWN"] + list(self.kwargs.get(
                "user_incar_settings", {}).keys()):
            # For these parameters as well as user specified settings, override
            # the incar settings.
            if parent_incar.get(k, None) is not None:
                incar[k] = parent_incar[k]
            else:
                incar.pop(k, None)

        # use new LDAUU when possible b/c the Poscar might have changed
        # representation
        if incar.get('LDAU'):
            u = incar.get('LDAUU', [])
            j = incar.get('LDAUJ', [])
            if sum([u[x] - j[x] for x, y in enumerate(u)]) > 0:
                for tag in ('LDAUU', 'LDAUL', 'LDAUJ'):
                    incar.update({tag: parent_incar[tag]})
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in incar:
                incar.update({"LMAXMIX": parent_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        incar["EDIFF"] = min(incar.get("EDIFF", 1), parent_incar["EDIFF"])
        return incar

    @property
    def kpoints(self):
        self._config_dict["KPOINTS"]["grid_density"] = self.grid_density
        kpoints = super(PRLStaticSet, self).kpoints
        # Prefer to use k-point scheme from previous run
        if self.prev_kpoints and self.prev_kpoints.style != kpoints.style:
            if self.prev_kpoints.style == Kpoints.supported_modes.Monkhorst:
                k_div = [kp + 1 if kp % 2 == 1 else kp
                         for kp in kpoints.kpts[0]]
                kpoints = Kpoints.monkhorst_automatic(k_div)
            else:
                kpoints = Kpoints.gamma_automatic(kpoints.kpts[0])
        return kpoints

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, grid_density=8000,
                       small_gap_multiply=None, **kwargs):
        """
        Generate a set of Vasp input files for static calculations from a
        directory of previous Vasp run.

        Args:
            prev_calc_dir (str): Directory containing the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            standardize (float): Whether to standardize to a primitive
                standard cell. Defaults to False.
            sym_prec (float): Tolerance for symmetry finding. If not 0,
                the final structure from the previous run will be symmetrized
                to get a primitive standard cell. Set to 0 if you don't want
                that.
            international_monoclinic (bool): Whether to use international
                    convention (vs Curtarolo) for monoclinic. Defaults True.
            grid_density (int): density of k-mesh by reciprocal atom (defaults to 8000)
            small_gap_multiply ([float, float]): If the gap is less than
                1st index, multiply the default reciprocal_density by the 2nd
                index.
            \\*\\*kwargs: All kwargs supported by MPStaticSet,
                other than prev_incar and prev_structure and prev_kpoints which
                are determined from the prev_calc_dir.
        """
        vasprun, outcar = get_vasprun_outcar(prev_calc_dir)

        prev_incar = vasprun.incar
        prev_kpoints = vasprun.kpoints

        # We will make a standard structure for the given symprec.
        prev_structure = get_structure_from_prev_run(
            vasprun, outcar, sym_prec=standardize and sym_prec,
            international_monoclinic=international_monoclinic)

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                grid_density = grid_density * small_gap_multiply[1]

        return PRLStaticSet(structure=prev_structure, prev_incar=prev_incar,
            prev_kpoints=prev_kpoints, grid_density=grid_density, **kwargs)
