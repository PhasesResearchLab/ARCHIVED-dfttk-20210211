"""
This module implements the Quasi-harmonic Debye approximation that can
be used to compute thermal properties.

Quasiharmonic Debye-Gruneisen model based on pymatgen's QHA and modified
to clean up and use the proper Gruneisen parameter.
"""
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

from collections import defaultdict

import numpy as np

from scipy.constants import physical_constants
from scipy.integrate import quadrature
from scipy.misc import derivative
from scipy.optimize import minimize

from pymatgen.core.units import FloatWithUnit
from pymatgen.analysis.eos import EOS, PolynomialEOS
from pymatgen.io.vasp.outputs import Vasprun

from prlworkflows.analysis.thermal_electronic import calculate_thermal_electronic_contribution

__author__ = "Kiran Mathew, Brandon Bocklund"
__credits__ = "Cormac Toher"


class QuasiharmonicDebyeApprox(object):
    """
    Args:
        energies (list): list of DFT energies in eV
        volumes (list): list of volumes in Ang^3
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        eos (str): equation of state used for fitting the energies and the
            volumes.
            options supported by pymatgen: "quadratic", "murnaghan", "birch",
                "birch_murnaghan", "pourier_tarantola", "vinet",
                "deltafactor", "numerical_eos"
        pressure (float): in GPa, optional.
        poisson (float): poisson ratio.
        anharmonic_contribution (bool): Whether to use the Debye-Gruneisen model
    """
    def __init__(self, energies, volumes, structure, vasprun_paths, t_min=0.0, t_step=10,
                 t_max=2000.0, eos="vinet", pressure=0.0, poisson=0.25,
                 gruneisen=True, bp2gru=2./3.):
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        self.temperature_min = t_min
        self.temperature_max = t_max
        self.temperature_step = t_step
        self.eos_name = eos
        self.pressure = pressure
        self.poisson = poisson
        self.bp2gru = bp2gru
        self.gruneisen = gruneisen
        self.mass = sum([e.atomic_mass for e in self.structure.species])
        self.natoms = self.structure.composition.num_atoms
        self.avg_mass = physical_constants["atomic mass constant"][0] * self.mass / self.natoms  # kg
        self.kb = physical_constants["Boltzmann constant in eV/K"][0]
        self.hbar = physical_constants["Planck constant over 2 pi in eV s"][0]
        self.gpa_to_ev_ang = 1./160.21766208  # 1 GPa in ev/Ang^3
        self.gibbs_free_energy = []  # optimized values, eV
        # list of temperatures for which the optimized values are available, K
        self.temperatures = []
        self.optimum_volumes = []  # in Ang^3
        # fit E and V and get the bulk modulus(used to compute the Debye
        # temperature)
        print("Fitting E and V")
        self.eos = EOS(eos)
        self.ev_eos_fit = self.eos.fit(volumes, energies)
        self.bulk_modulus = self.ev_eos_fit.b0_GPa  # in GPa

        vasprun_objs = [Vasprun(v) for v in vasprun_paths]
        thermal_electronic_props = [calculate_thermal_electronic_contribution(vr.complete_dos, t0=t_min, t1=t_max, td=t_step, natom=1) for vr in vasprun_objs]
        self.F_el = [p['free_energy'] for p in thermal_electronic_props]
        self.F_vib = [np.zeros(self.F_el[0].size) for _ in volumes]
        self.optimize_gibbs_free_energy()

    def optimize_gibbs_free_energy(self, verbose=False):
        """
        Evaluate the gibbs free energy as a function of V, T and P i.e
        G(V, T, P), minimize G(V, T, P) wrt V for each T and store the
        optimum values.

        Note: The data points for which the equation of state fitting fails
            are skipped.
        """
        temperatures = np.linspace(
            self.temperature_min,  self.temperature_max,
            int(np.ceil((self.temperature_max - self.temperature_min)
            / self.temperature_step) + 1))

        for t in temperatures:
            G_opt, V_opt = self.optimizer(t, verbose=verbose)
            try:
                #G_opt, V_opt = self.optimizer(t)
                pass
            except:
                if len(temperatures) > 1:
                    print("EOS fitting failed, so skipping this data point, {}".
                          format(t))
                    continue
                else:
                    raise
            self.gibbs_free_energy.append(G_opt)
            self.temperatures.append(t)
            self.optimum_volumes.append(V_opt)

    def optimizer(self, temperature, verbose=False):
        """
        Evaluate G(V, T, P) at the given temperature(and pressure) and
        minimize it wrt V.

        1. Compute the  vibrational helmholtz free energy, A_vib.
        2. Compute the gibbs free energy as a function of volume, temperature
            and pressure, G(V,T,P).
        3. Preform an equation of state fit to get the functional form of
            gibbs free energy:G(V, T, P).
        4. Finally G(V, P, T) is minimized with respect to V.

        Args:
            temperature (float): temperature in K

        Returns:
            float, float: G_opt(V_opt, T, P) in eV and V_opt in Ang^3.
        """
        G_V = []  # G for each volume
        # G = E(V) + PV + A_vib(V, T) + F_el(V, T)
        for i, v in enumerate(self.volumes):
            vib = self.vibrational_free_energy(temperature, v)
            temp_idx = int((temperature-self.temperature_min)/self.temperature_step)
            self.F_vib[i][temp_idx] = vib
            e_wo_el = self.energies[i] + self.pressure * v * self.gpa_to_ev_ang + vib
            el = self.F_el[i][temp_idx]
            G_V.append( e_wo_el + el)
            #G_V.append(e_wo_el)
            if verbose:
                print((temperature - self.temperature_min) / self.temperature_step)
                print('energy wo electronic {}'.format(e_wo_el))
                print('vibrational energy {}'.format(vib))
                print('energy electronic {}'.format(el))

        # fit equation of state, G(V, T, P)
        eos_fit = self.eos.fit(self.volumes, G_V)
        # minimize the fit eos wrt volume
        # Note: the ref energy and the ref volume(E0 and V0) not necessarily
        # the same as minimum energy and min volume.
        volume_guess = eos_fit.volumes[np.argmin(eos_fit.energies)]
        min_wrt_vol = minimize(eos_fit.func, volume_guess)
        # G_opt=G(V_opt, T, P), V_opt
        return min_wrt_vol.fun, min_wrt_vol.x[0]

    def vibrational_free_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9./8. * y + 3 * np.log(1 - np.exp(-y)) - self.debye_integral(y))

    def vibrational_internal_energy(self, temperature, volume):
        """
        Vibrational internal energy, U_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: vibrational internal energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9./8. * y +
                                                      3*self.debye_integral(y))

    def debye_temperature(self, volume):
        """
        Calculates the debye temperature.
        Eq(6) in doi.org/10.1016/j.comphy.2003.12.001. Thanks to Joey.

        Eq(6) above is equivalent to Eq(3) in doi.org/10.1103/PhysRevB.37.790
        which does not consider anharmonic effects. Eq(20) in the same paper
        and Eq(18) in doi.org/10.1016/j.commatsci.2009.12.006 both consider
        anharmonic contributions to the Debye temperature through the Gruneisen
        parameter at 0K (Gruneisen constant).

        The anharmonic contribution is toggled by setting the anharmonic_contribution
        to True or False in the QuasiharmonicDebyeApprox constructor.

        Args:
            volume (float): in Ang^3

        Returns:
            float: debye temperature in K
         """
        term1 = (2./3. * (1. + self.poisson) / (1. - 2. * self.poisson))**1.5
        term2 = (1./3. * (1. + self.poisson) / (1. - self.poisson))**1.5
        f = (3. / (2. * term1 + term2))**(1. / 3.)
        debye = 2.9772e-11 * (volume / self.natoms) ** (-1. / 6.) * f * np.sqrt(self.bulk_modulus/self.avg_mass)
        if self.gruneisen:
            # bp2gru should be the correction to the Gruneisen constant.
            # High temperature limit: 2/3
            # Low temperature limit: 1
            # take 0 K E-V curve properties
            dBdP = self.ev_eos_fit.b1  # bulk modulus/pressure derivative
            gamma = (1+dBdP)/2 - self.bp2gru  # 0K equilibrium Gruneisen parameter
            return debye * (self.ev_eos_fit.v0 / volume) ** (gamma)
        else:
            return debye


    @staticmethod
    def debye_integral(y):
        """
        Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

        Args:
            y (float): debye temperature/T, upper limit

        Returns:
            float: unitless
        """
        # floating point limit is reached around y=155, so values beyond that
        # are set to the limiting value(T-->0, y --> \infty) of
        # 6.4939394 (from wolfram alpha).
        factor = 3. / y ** 3
        if y < 155:
            integral = quadrature(lambda x: x ** 3 / (np.exp(x) - 1.), 0, y)
            return list(integral)[0] * factor
        else:
            return 6.493939 * factor


    def get_summary_dict(self):
        """
        Returns a dict with a summary of the computed properties.
        """
        d = defaultdict(list)
        d["pressure"] = self.pressure
        d["poisson"] = self.poisson
        d["mass"] = self.mass
        d["natoms"] = int(self.natoms)
        d["bulk_modulus"] = self.bulk_modulus
        d["gibbs_free_energy"] = self.gibbs_free_energy
        d["temperatures"] = self.temperatures
        d["optimum_volumes"] = self.optimum_volumes
        return d
