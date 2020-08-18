# -------- energy in eV, temperature in K
# assume every variable starting with a-h  and o-z are real numbers
# common block named comcon
from __future__ import division
import sys
import os
import subprocess
import math
import copy
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import brentq, curve_fit
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
from scipy.interpolate import UnivariateSpline
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import dfttk.pyphon as ywpyphon
from dfttk.utils import sort_x_by_y


k_B = physical_constants['Boltzmann constant in eV/K'][0]

def substr(str1, str2, pos):
  try:
    if str1.index(str2)==pos:
        #print("idx=",str1.index(str2))
        return True
    else:
        return False
  except ValueError:
    return False

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False


def VBMinterpolate(v,e,dos,mode='linear'):
  # obselete, not used
  x = np.zeros(9)
  y = np.zeros(9)
  x[0:8]=e
  x[8] = 0.0
  y[0:8]=dos
  y[8] = 0.0
  f = interp1d(x,y,kind=mode)
  print(v,f(v))
  return(f(v))
  
# this is a FORTRAN function (e.g. 1 return value)
def pregetdos(f): # Line 186
    """
    to make the code can also handle WIEN2k dos in the unit of eV

    Parameters
    ----------
    f : file descriptor for the DOS file

    Returns
    -------
    xdn : lower energy to integrate over?
    xup : higher energy to integrate over?
    vde : band energy intercal
    e (array): band energy mesh the Fermi energy has been shifted to zero
    DOS (array) : e dos
    """
    # read the file
    lines = f.readlines() # read in all lines then determine is it is WIEN2k DOS (in the unit eV) file or VASP DOS file
    # now the first line should be the one with the data, lets remove the one into its own special line
    tmp = lines[0]
    if substr(tmp,"#  BAND", 0):  
        tmp = lines[1]
        tmp1 = lines[2]
        if substr(tmp, "#EF=",0) and substr(tmp1, "# ENERGY",0):  
            tmp1 = tmp[31:43].replace("NENRG=","")
            if isint(tmp1):
                n_dos = int(tmp1)
                tmp = lines[2]
                lines = lines[3:n_dos+3]
                wienEdos = np.zeros(n_dos)
                ve = np.zeros(n_dos)
                for i, l in enumerate(lines):
                    split_l = l.split(' ')
                    split_l = [k for k in split_l if k != '']
                    ve[i], wienEdos[i] = (float(split_l[0]), float(split_l[1]))
                edn = ve[0]
                eup = ve[n_dos-1]
                ve = np.linspace(edn, eup, n_dos)
                vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?
                return edn, eup, vde, ve, wienEdos

    tmp = lines[5]
    data_line = tmp[0:32].split(' ') #n_dos >10000, no space left before it in VASP
    data_line.extend(tmp[32:].split(' '))
    # filter out empty spaces
    data_line = [k for k in data_line if k != '']
    #print (data_line)
    eup, edn, n_dos, eFermi = (float(data_line[0]),
                           float(data_line[1]),
                           int(data_line[2]),
                           float(data_line[3])) # we're leaving the last number behind
    lines = lines[6:n_dos+6]
    # line 197 goes to line 209
    eup = eup - eFermi
    edn = edn - eFermi
    vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?

    # vectors
    ve = np.linspace(edn, eup, n_dos)
    vaspEdos = np.zeros(n_dos)

    for i, l in enumerate(lines):
        # why do we need to do this?
        split_l = l.split(' ')
        # filter again
        split_l = [k for k in split_l if k != '']
        if len(split_l)>=5: #spin polarized
            t, vaspEdos[i], y, vdos, x = (float(split_l[0]), float(split_l[1]), float(split_l[2]), float(split_l[3]), float(split_l[4]))
            vaspEdos[i] += y
        else:
            t, vaspEdos[i], vdos = (float(split_l[0]), float(split_l[1]), float(split_l[2]))
    return edn, eup, vde, ve, vaspEdos

def getdos(xdn, xup, dope, NEDOS, gaussian, edn, eup, vde, ve, tdos): # Line 186
    """

    Parameters
    ----------
    dos : pymatgen.electronic_structure.dos.Dos
        DOS object from pymatgen
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    dos_grid_size : int
        Number of DOS points have in the energy/density grid
    gaussian_grid_size : int
        Number of Gaussian points to use in the grid mesh around the Fermi energy

    Returns
    -------
    tuple
        Tuple of a (float, float, array, array) of the number of electrons,
        Fermi level shift due to doping, and the arrays of energies and densities on a grid.
    """

    for i,energy in enumerate(ve):
      if energy <-15.0 and tdos[i]==0.0:
        xdn = energy

    n_dos = len(tdos)
    idx = closest(ve,0.0)
    for i in range(idx,n_dos):
        if tdos[i]!=0.0:
            iBoF = i-1
            break

    eBoF = -1.0
    if iBoF>=idx:
      eBoF = ve[iBoF]
      espr = tdos[iBoF+2]-tdos[iBoF+1]
      if espr>0.0:
        espr = tdos[iBoF+1]/espr*vde
        if (espr < vde):
          eBoF = ve[iBoF+1] - espr

    #print("eBoF=", eBoF)
    xdn = max(xdn,edn)
    xup = min(xup,eup)

    e = np.linspace(xdn,xup,NEDOS,dtype=float)
    if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, 0.0, eBoF, NEDOS)

    dos = refdos(eBoF, 0.0, vde, edn, e, ve, tdos)
    ados = cumtrapz(dos, e, initial=0.0)
    idx = closest(e,0.0)
    for idx1 in range(idx-1, 0, -1):
        if ados[idx1] != ados[idx] : break
    NELECTRONS = ados[idx] - e[idx]/(e[idx1]-e[idx])*(ados[idx1]-ados[idx])

    dF = 0.0
    if dope != 0.0:
        NELECTRONS = NELECTRONS+dope
        idx = closest(ados,NELECTRONS)
        for idx1 in range(idx-1, 0, -1):
            if ados[idx1] != ados[idx] : break
        #if idx == (NEDOS-1) or ados[idx] == ados[idx+1]:
        #print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
        if idx1 <= 0 or idx >= (NEDOS-1) or ados[idx] == ados[idx1]:
            print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
            # we are dopidxng too much
            raise ValueError('Too much doping')
        dF = (NELECTRONS-ados[idx])/(ados[idx1] - ados[idx])*(e[idx1] - e[idx])+e[idx]
                # dF is the shift in the Fermi energy due to doping
        e = e - dF # This is done in a loop (line 289), but I think we can do without

    if gaussian != 0.0 and abs(dope)>0.0001: # why did I do this ***********************
    #if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, dF, eBoF, NEDOS)

    dos = refdos(eBoF, dF, vde, edn, e, ve, tdos)
    edos = e*dos
    ados = cumtrapz(dos, e, initial=0.0)
    energy = cumtrapz(edos, e, initial=0.0)
    idx = closest(e,0.0)
    NELECTRONS = ados[idx] - e[idx]/(e[idx+1]-e[idx])*(ados[idx+1]-ados[idx])
    E0 = energy[idx] - e[idx]/(e[idx+1]-e[idx])*(energy[idx+1]-energy[idx])

    return NELECTRONS, E0, dF, e, dos, eBoF

def remesh(xdn, xup, gaussian, dF, eBoF, NEDOS):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameters
    ----------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    ve : original e mesh
    NEDOS : original e mesh
    gaussian : parameter used to refine the e mesh near the Fermi energy

    Return
    ------
    e : refined e mesh
    """

    e = np.zeros(NEDOS)
    e[0] = xdn - dF
    xde = 2.0*(xup - xdn)/(NEDOS-1)
    if eBoF>0.0:
        xde = 3.0*(xup - xdn)/(NEDOS-1)
    sigma = -0.5*(gaussian/(xup-xdn))**2
    fac = gaussian/(math.sqrt(2.0*math.pi))
    for i in range(1,NEDOS):
        f1 = 1.0 + fac*math.exp(sigma*(e[i-1])**2)
        if eBoF>0.0:
          if dF < eBoF:
             f1 += fac*math.exp(sigma*((e[i-1]-eBoF+dF))**2)
          else:
             f1 += fac*math.exp(sigma*((e[i-1]+dF))**2)
        e[i] = e[i-1]+xde/f1
    return e

def refdos(eBoF, dF, vde, edn, e, ve, tdos):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameter
    ---------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    e : refined e mesh
    ve : original e mesh
    tdos : original e dos
    Return
    ------
    dos : refined e dos
    """

    dos = np.zeros(len(e))
    n_dos = len(tdos)
    for i in range(0, len(e)):
        tx = e[i] + dF
        kx = int((tx-edn)/vde) # Converted to int, remember the type!
        kx = max([kx,0]) # XXX: is this translated correctly? What is the 1 in fortran?
        kx = min([n_dos-2, kx]) # TODO: the ndos-1 was here before. could be a source of error
        if tdos[kx+1]==0.0 and ve[kx+1]>0.0 and ve[kx+1]<vde:
          # handling near the Top of valence band
          if tx >= 0.0:
            dos[i] = 0.0
          else:
            dos[i] = tdos[kx]*tx/ve[kx]
            #dos[i] = tdos[kx]*(tx/ve[kx])**2
            #dos[i] = VBMinterpolate(tx,ve[kx-8:kx],tdos[kx-8:kx],mode='quadratic')
        elif eBoF > 0.0 and tdos[kx]==0.0 and ve[kx+1]-eBoF<vde and ve[kx+1]-eBoF>0.0:
          # handling near the bottom of conduction band
            if tx <= eBoF:
              dos[i] = 0.0
            else:
              dos[i] = tdos[kx+1]*(tx-eBoF)/(ve[kx+1]-eBoF)
        else:
          dos[i] = tdos[kx] + (tdos[kx+1] - tdos[kx])/vde*(tx - ve[kx])
    return dos
    
def closest(e,val):
    """
    find the index of the band energy which is the close to the energy val
    Parameters
    ----------
    e : float
    array of band energy for the e dos
    val : given value of band energy

    Return
    ------
    index of e that closest to the energy val
    """
    idx = np.abs(e-val).argmin()
    if e[idx] < val:
        idx = idx + 1
    return idx

def gfind(mu_el, pe, pdos, NELECTRONS, Beta, IntegrationFunc=trapz):
    """
    Calculate the number of electron difference from 0K given chemical potential. the purpose is the find the
    chemical potential to make zero of number of electron difference from 0K

    Parameters
    ----------
    mu_el : chemical potential, :math:`\mu`, in the Fermi distribution
    pe : eigenenergies
    pdos : density of states (:math:`n(\varepsilon) \varepsilon`
    NELECTRONS : Total number of electrons in the system at 0K
    Beta : :math:`\frac{1}{T*k_{B}}`

    Returns
    -------
    The number of electron difference from 0K given chemical potential
    """

    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k = len(tc)
    fn = pdos[0:k]/(np.exp(tc[0:k])+1.0)
    return IntegrationFunc(fn, pe[0:k])- NELECTRONS


# line 363
def caclf(pe, pdos, NELECTRONS, Beta, mu_ref=0.0, dF=0.0, IntegrationFunc=trapz): #line 363
    """
    Calculate thermal free energy from electronic density of states (e DOS)

    Parameters
    ----------
    pe : band energy array
    pdos : e DOS
    NELECTRONS : total number of electrons
    Beta : 1/(kB*T)

    Returns
    -------
    electron chememical potential, internal energy, entropy, carrier amount, coefficient to cal Seebeck
    """

    #print ("dF=", dF)
    if 1==1:
        mu_el = brentq(gfind, mu_ref-10.0, mu_ref+10.0, args=(pe, pdos, NELECTRONS, Beta, IntegrationFunc), maxiter=10000)
    else:
        t0 = mu_ref
        d0 = gfind(t0, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
        if d0 > 0.0: td = -0.1
        elif d0 <0.0: td = 0.1
        else: return t0
        for i in range(999):
            t1 = t0 + td
            d1 = gfind(t1, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d1*d0 < 0.0: break
            elif d1*d0 == 0.0: break
            t0 = t1
            d0 = d1
            td = td + td
        for i in range(999):
            t2 = (t0 + t1)*0.5
            d2 = gfind(t2, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d2*d0 < 0.0:
                t1 = t2
                d1 = d2
            else:
                t0 = t2
                d0 = d2
            if abs(t1-t0) <1.e-8:
                mu_el = 0.5*(t0+t1)
                break
        #mu_el = brentq(gfind, t0, t1, args=(pe, pdos, NELECTRONS, Beta), maxiter=10000)
        print("xxxxxxx", mu_el,mu_old)
    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k1 = len(tc)
    tf = 1.0/(np.exp(tc)+1.0)
    fn = pdos[0:k1]*pe[0:k1]*tf
    u = IntegrationFunc(fn, pe[0:k1])

    k0 = closest(tc,-200)
    tf0 = tf[k0:]
    pdos = pdos[k0:k1]
    pe = pe[k0:k1]
    tf1 = 1.0 - tf0 + 1.e-60 # 1.e-60 is used to avoid log exception
    fn = pdos*(tf0*np.log(tf0)+tf1*np.log(tf1))
    s = IntegrationFunc(fn, pe)

    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    fn2 = pdos*tf*(pe-mu_el)
    Q_el = IntegrationFunc(fn, pe)
    Q_p = IntegrationFunc(fn[pe<=dF], pe[pe<=dF])
    Q_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])
    Y_el = IntegrationFunc(fn2, pe)

    fn = pdos*(pe-mu_el)*tf
    if Q_el!=0.0:
        e_ = IntegrationFunc(fn, pe)/Q_el
        fn = pdos[0:k1]*(pe[0:k1]-mu_el-e_)**2*tf
        cv = IntegrationFunc(fn, pe[0:k1])
    else:
        cv = 0.0
    fn = pdos[0:k1]*(pe[0:k1]-mu_el)**2*tf
    c_mu = IntegrationFunc(fn, pe[0:k1])

#   hole/electron concentration by effective carrier
    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    f2 = interp1d(pe, fn, kind='linear')
    fmu = f2(mu_el)
    x = np.hstack([pe[pe<mu_el],mu_el])
    y = np.hstack([fn[pe<mu_el],fmu])
    W_p = IntegrationFunc(y,x)
    x = np.hstack([mu_el, pe[pe>mu_el]])
    y = np.hstack([fmu, fn[pe>mu_el]])
    W_e = IntegrationFunc(y,x)
    #W_e = IntegrationFunc(fn[pe>mu_el], pe[pe>mu_el])
    #W_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])

#   hole/electron concentration by alternative difination
    fn = pdos*(1.0-tf0)
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([pe[pe<mu_el],mu_el])
    #y = np.hstack([fn[pe<mu_el],fmu])
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([pe[pe<dF],dF])
    y = np.hstack([fn[pe<dF],fmu])
    Y_p = IntegrationFunc(y,x)

    fn = pdos*tf0
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([mu_el, pe[pe>mu_el]])
    #y = np.hstack([fmu, fn[pe>mu_el]])
    #print ("mu_el", mu_el, dF)
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([dF, pe[pe>dF]])
    y = np.hstack([fmu, fn[pe>dF]])
    Y_e = IntegrationFunc(y,x)

    return mu_el, u, -s*k_B, cv*k_B*Beta*Beta, Q_el, Y_el, Q_p, Q_e, c_mu*k_B*Beta*Beta, W_p, W_e, Y_p, Y_e


def T_remesh(t0, t1, td, _nT=-1):
    T = []
    if td > 0:
        for t in np.arange(t0,t1+td, td):
            T.append(t)
        return np.array(T)

    if _nT <= 0: nT = 51
    else: nT = _nT
    a = 100./nT
    dT_new = abs(td)/(1+(nT-1)*0.5*a)
    for i in range (nT):
      T.append(t0+i*dT_new*(1+i*a))
    T = np.array(T)
    p = (t1-t0)/(max(T)-t0)
    for i in range (nT):
      T[i] = round((T[i]-t0)*p+t0,2)
    return T


def runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, 
    _T=[], dos=sys.stdin, fout=sys.stdout, vol=None, IntegrationFunc=trapz):
    """
    Calculate thermal free energy from electronic density of states (e DOS)

    Parameters
    ----------
    t0 : float
        Low temperature limit
    t1 : float
        High temperature limit
    td : float
        Temperature increment
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize
        the calculated properties in the unit of per atom
    dos : file description for the DOSCAR or pymatgen dos object
        Filename for VASP DOSCAR
    outf : file description
        Output file description for the calculated properties

    Return
    ------
    Tuple of 14 float array containing
    thermal electron free energy, entropy, specific heat, M_el, seebeck_coefficients, 
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, temperature.
    Other quantities are for researching purpose
    """

    if hasattr(dos, 'read'):
        edn, eup, vde, dos_energies, vaspEdos = pregetdos(dos) # Line 186
    else:
        e_fermi = dos.efermi
        eup = np.max(dos.energies) - e_fermi
        edn = np.min(dos.energies) - e_fermi
        n_dos = len(dos.energies) # number of points in DOS
        vde = (eup - edn)/(n_dos-1) # change in energy per step
        dos_energies = np.linspace(edn, eup, n_dos) # linearize: sometimes rounding errors in DOSCAR
        vaspEdos = np.array(dos.get_densities())
    NELECTRONS, E0, dF, e, dos, Eg = getdos(xdn, xup, dope, ndosmx, gaussian, edn, eup, vde, dos_energies, vaspEdos)

    if Eg < 0.0: Eg = 0.0
    if vol == None:
        fout.write('#Bandgap= {} eV. '.format(Eg))
    else:
        fout.write('#Bandgap= {} eV at volume= {} Angstrom^3/cell. '.format(Eg,vol))

    fout.write('Fermi energy was shifted {} due to doping of {} resulting Ne={} \n'.format(dF, dope, NELECTRONS))

    # for all temperatures
    if len(_T)!=0:
      T=copy.deepcopy(_T)
    elif td>0.0:
      T = np.arange(t0,t1+td,td) # temperature
    else:
      if self.debug:
        T = T_remesh(t0,t1,td,_nT=65)
      else:
        T = T_remesh(t0,t1,td,_nT=129)
    nT = len(T)
    U_el = np.zeros(nT)
    S_el = np.zeros(nT)
    C_el = np.zeros(nT) # electronic specific heat
    C_mu = np.zeros(nT) # electronic specific heat at constant chemical potential
    M_el = np.zeros(nT) # electronic chemical potential, i.e., absolute thermal electric force
    Q_el = np.zeros(nT) # total number of thermal Carrier
    Y_el = np.zeros(nT)
    Q_p = np.zeros(nT)
    Q_e = np.zeros(nT)
    W_p = np.zeros(nT)
    W_e = np.zeros(nT)
    Y_p = np.zeros(nT)
    Y_e = np.zeros(nT)
    seebeck_coefficients = np.zeros(nT)
    U_el[0] = E0

    for i in range(0,nT):
        if T[i]==0.0: continue
        Beta = 1.0e0/(T[i]*k_B)
        M_el[i], U_el[i], S_el[i], C_el[i], Q_el[i],Y_el[i], Q_p[i],Q_e[i],  C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i] = caclf(e, dos, NELECTRONS, Beta, M_el[i-1], dF=-dF, IntegrationFunc=IntegrationFunc)
        if Q_el[i]>0.0:
            seebeck_coefficients[i] = -1.0e6*Y_el[i]/Q_el[i]/T[i]

    F_el_atom = (U_el - T * S_el - E0) / natom  # electronic free energy per atom
    S_el_atom = S_el / natom  # entropy per atom
    #dU_dT = np.gradient(U_el, td) # gradient on U_el with step size of td
    #C_el_atom = dU_dT/natom # electronic heat capacity per atom
    C_el_atom = C_el/natom # electronic heat capacity per atom

    return F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e


def thelecAPI(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, doscar):
    """
    API to calculate the thermal electronic properties from DOSCAR

    Parameters
    ----------
    t0 : float
        Low temperature limit
    t1 : float
        High temperature limit
    td : float
        Temperature increment
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize
        the calculated properties in the unit of per atom
    doscar : str
        Filename for VASP DOSCAR

    Output to (printed to outf)
    ---------------------------
    The properties in the order of temperature, thermal electron free energy, entropy, 
    specific heat, seebeck_coefficients, Lorenz number,  
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat
    """

    with open(doscar, 'r') as fp:
      with open(outf, 'w') as fvib:
        F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=fp, fout=fvib)
        fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
        for i in range(T.size):
            L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
            if Q_el[i] != 0.0: L = C_el_atom[i]/Q_el[i]*k_B
            if Q_el[i] > 1.e-16: L = C_el_atom[i]/Q_el[i]*k_B
            fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i]))

def BMvol(V,a):
  T = V**(-1./3)
  fval = a[0]+a[1]*T
  if len(a) > 2:
    fval += a[2]*T*T
  if len(a) > 3:
    fval += a[3]*T*T*T
  if len(a) > 4:
    fval += a[4]*T*T*T*T
  return(fval)

def BMvolP(V,a):
  T = V**(-1./3)
  fval = a[1]
  if len(a) > 2:
    fval += 2*a[2]*T
  if len(a) > 3:
    fval += 3*a[3]*T*T
  if len(a) > 4:
    fval += 4*a[4]*T*T*T
  return(fval)

def BMvolB(V,a):
  T = V**(-1./3)
  fval = 2*a[2]
  if len(a) > 3:
    fval += 2*3*a[3]*T
  if len(a) > 4:
    fval += 3*4*a[4]*T*T
  return(fval)

def BMvol4(T, a, b, c, d):
  return (BMvol(T, [a,b,c,d]))

def BMvol5(T, a, b, c, d, e):
  return (BMvol(T, [a,b,c,d,e]))

def BMfitB(V, x, y, BMfunc):
  f, pcov = curve_fit(BMfunc, x, y)
  p = BMvolP(V, f)
  b = BMvolB(V, f)
  P = p*V**(-4./3)*(-1./3.)
  B = b*V**(-4./3)*(-1./3.)*V**(-4./3)*(-1./3.) + p*V**(-7./3)*(-1./3.)*(-4./3.)
  return B*V, P

def BMfitP(V, x, y, BMfunc):
  f, pcov = curve_fit(BMfunc, x, y)
  p = BMvolP(V, f)
  b = BMvolB(V, f)
  P = p*V**(-4./3)*(-1./3.)
  return P

def BMfitF(V, x, y, BMfunc):
  f, pcov = curve_fit(BMfunc, x, y)
  return BMvol(V,f)

def BMsmooth(_V, _E0, _Flat, _Fel, _Slat, _Sel, BMfunc):
    E0 = BMfitF(_V, _V, _E0, BMfunc)
    Flat = UnivariateSpline(_V, _Flat)(_V)
    Fel = UnivariateSpline(_V, _Fel)(_V)
    Slat = UnivariateSpline(_V, _Slat)(_V)
    Sel = UnivariateSpline(_V, _Sel)(_V)
    return E0, Flat, Fel, Slat, Sel

def CenDif(v, vol, F, N=7,kind='cubic'):
    dV = 0.01*(max(vol) - min(vol))
    if kind=='cubic':
        f2 = interp1d(vol, F, kind='cubic')
    else:
        f2 = UnivariateSpline(vol, F)
    result = 0.0 
    n = N//2
    for i in range(0, n):
        #print (v, f2(v+dV*(i+1)) , f2(v-dV*(i+1)))
        result += (f2(v+dV*(i+1)) - f2(v-dV*(i+1)))/(2.0*(i+1)*dV)
    #sys.exit()
    return result/n


def CenDifB(_v, vol, F, BMfunc, N=7,kind='cubic'):
    dV = 0.01*(max(vol) - min(vol))
    try: 
    #if True:
      if _v!=0.0:
          v = brentq(BMfitP, _v-N*dV, _v+N*dV, args=(vol, F, BMfunc), maxiter=10000)
      else:
          v = brentq(BMfitP, min(vol), max(vol), args=(vol, F, BMfunc), maxiter=10000)
      #ff = BMfitF(v, vol, F, BMfunc)
      ff = interp1d(vol, F)(v)
    except:
      v = _v
      ff = None
    result = 0.0 
    n = N//2
    for i in range(0, n):
        result += (CenDif(v+dV*(i+1), vol, F, kind=kind) - CenDif(v-dV*(i+1), vol, F, kind=kind))/(2.0*(i+1)*dV)
    return result*v/(n), v, ff

def BMDifB(_v, vol, F, BMfunc, N=7):
    dV = 0.01*(max(vol) - min(vol))
    if _v!=0.0:
        v = brentq(BMfitP, _v-N*dV, _v+N*dV, args=(vol, F, BMfunc), maxiter=10000)
    else:
        v = brentq(BMfitP, min(vol), max(vol), args=(vol, F, BMfunc), maxiter=10000)

    #ff = BMfitF(v, vol, F, BMfunc)
    ff = interp1d(vol, F)(v)
    bb, pp = BMfitB(v, vol, F, BMfunc)
    return bb, v, ff, pp


def SymDif(v, vol, F, N=7,kind='cubic'):
    dV = 0.001*(max(vol) - min(vol))
    if kind=='cubic':
        f2 = interp1d(vol, F, kind='cubic')
    else:
        f2 = UnivariateSpline(vol, F)
    result = 0.0
    for i in range(0, N+1):
        result += (f2(v+dV*(i)) - f2(v-dV*(i+N)))/(N*dV)
    return result/(N+1)


def SymDifB(v, vol, F, N=7,kind='cubic'):
    dV = 0.001*(max(vol) - min(vol))
    result = 0.0
    for i in range(0, N+1):
        result += (SymDif(v+dV*(i), vol, F, kind=kind) - SymDif(v-dV*(i+N), vol, F, kind=kind))/(N*dV)
    return result*v/(N+1)


def debye_heat_capacity(temperature, debye_T, natoms):
    """
    debye Vibrational heat capacity, C_vib(V, T).
    Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

    Args:
        temperature (float): temperature in K
        volume (float)

    Returns:
        float: vibrational heat capacity in eV
    """
    y = debye_T / temperature
    factor = 3. / y ** 3
    if y < 155:
        integral = quadrature(lambda x: x ** 4 *np.exp(x)/ (np.exp(x) - 1.)**2, 0, y)
        return 3*k_B * natoms * list(integral)[0] * factor
    else:
        return k_B * natoms * 4./5.*math.pi**4 * factor

def debye_phonon(x, temperature, natoms, clat):
    return debye_heat_capacity(temperature, x, natoms) - clat

def get_debye_T_from_phonon_Cv(temperature, clat, dlat, natoms, _td=50):
    if temperature <=0: return dlat
    t0 = dlat
    d0 = debye_phonon(t0, temperature, natoms, clat)
    if d0 > 0.0: td = _td
    elif d0 <0.0: td = -_td 
    else: return t0
    for i in range(999):
        if t0 < 0.1 : return 0
        t1 = t0 + td
        t1 = max(t1,0.01)
        d1 = debye_phonon(t1, temperature, natoms, clat)
        if d1*d0 < 0.0: break
        t0 = t1
        d0 = d1
        td = td + td
    return brentq(debye_phonon, t0, t1, args=(temperature, natoms, clat), maxiter=10000)
        
    """
    for i in range(20):
        if abs(t1-t0) < 0.0001: break
        t2 = 0.5*(t0 + t1)
        d2 = debye_phonon(t2, temperature, natoms, clat)
        if d2*d0 < 0.0:
            t1 = t2
            d1 = d2
        elif d2*d0 > 0.0:
            t0 = t2
            d0 = d2
        else:
            return t2

    return 0.5*(t0 + t1)
    """

def vol_within(vol, volumes, thr=0.001):
    for i,v in enumerate(volumes):
        if (abs(vol-v) < thr*vol): return True
    return False


class thelecMDB():
    """
    API to calculate the thermal electronic properties from the saved dos and volume dependence in MongDB database

    Parameters
    ----------
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize 
        the calculated properties in the unit of per atom
    everyT : int
        Default 1. number of temperature points skipped from QHA analysis for debug speeding purpose 
    outf : str
        Output filename for the calculated properties
    db_file : str
        Filename containing the information to access the MongoDB database
    metatag : str
        metadata tag to access the to be calculated the compound
    qhamode : str
        Mode for the quasiharmonic calculations, according to it the T dependence volume to be extracted
    eqamode : int
        Mode to get LTC and the equilibrium volume

    Output to (printed to outf)
    ---------------------------
    The properties in the order of temperature, thermal electron free energy, entropy, 
    specific heat, seebeck_coefficients, Lorenz number,  
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, ..., 
    and the equilibrium volume extracted from MongoDB in the last column
    """

    def __init__(self, t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, db_file, 
                noel=False, everyT=1, metatag=None, qhamode=None, eqmode=0, elmode=1, smooth=False, plot=False,
                phasename=None, pyphon=False, debug=False, renew=False):
        from atomate.vasp.database import VaspCalcDb
        from pymatgen import Structure
        self.vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        self.t0 = t0
        self.t1 = t1
        self.td = td
        self.xdn = xdn
        self.xup = xup
        self.dope = dope
        self.ndosmx = ndosmx
        self.gaussian = gaussian
        self.natom = natom
        self.outf = outf
        self.noel = noel
        self.everyT = everyT
        self.tag = metatag
        self.qhamode = qhamode
        self.eqmode = eqmode
        self.elmode = elmode
        self.smooth = smooth
        self.plot = plot
        self.phasename=phasename
        self.pyphon=pyphon
        self.debug=debug
        self.renew=renew
        #print ("iiiii=",len(self._Yphon))


    def toYphon(self):
        self.Vlat = []
        self.Flat = []
        self.Clat = []
        self.Slat = []
        if self.debug:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=65)
        else:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=129)

        print ("extract the superfij.out used by Yphon ...")
        phdir = self.phasename+'/Yphon'
        if not os.path.exists(phdir):
            os.mkdir(phdir)
        for i in (self.vasp_db).db['phonon'].find({'metadata.tag': self.tag}):
            if vol_within(float(i['volume']), self.Vlat): continue
            self.Vlat.append(float(i['volume']))

            structure = Structure.from_dict(i['unitcell'])
            poscar = structure.to(fmt="poscar")
            unitcell_l = str(poscar).split('\n')
            natom = len(structure.sites)

            supercell_matrix = i['supercell_matrix']
            supercell_structure = copy.deepcopy(structure)
            supercell_structure.make_supercell(supercell_matrix)

            sa = SpacegroupAnalyzer(supercell_structure)
            primitive_unitcell_structure = sa.find_primitive()
            poscar = primitive_unitcell_structure.to(fmt="poscar")
            punitcell_l = str(poscar).split('\n')

            natoms = len(supercell_structure.sites)
            poscar = supercell_structure.to(fmt="poscar")
            supercell_l = str(poscar).split('\n')
            vol = 'V{:010.6f}'.format(float(i['volume']))
            voldir = phdir+'/'+vol
            if not os.path.exists(voldir):
               os.mkdir(voldir)
            structure.to(filename=voldir+'/POSCAR')
            with open (voldir+'/metadata.json','w') as out:
                mm = i['metadata']
                mm['volume'] = i['volume']
                out.write('{}\n'.format(mm))
            with open (voldir+'/superfij.out','w') as out:
                for line in range (2,5):
                    out.write('{}\n'.format(unitcell_l[line]))
                for line in range (2,5):
                    out.write('{}\n'.format(supercell_l[line]))
                out.write('{} {}\n'.format(natoms, natoms//natom))
                for line in range (7,natoms+8):
                    out.write('{}\n'.format(supercell_l[line]))
                force_constant_matrix = np.array(i['force_constants'])
                hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)
                for ii in range(natoms):
                    for jj in range(natoms):
                        for x in range(3):
                            for y in range(3):
                                hessian_matrix[ii*3+x, jj*3+y] = -force_constant_matrix[ii,jj,x,y]
                for xx in range(natoms*3):
                    for yy in range(natoms*3-1):
                        out.write('{} '.format(hessian_matrix[xx,yy]))
                    out.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))

            cwd = os.getcwd()
            os.chdir( voldir )
            if not os.path.exists('vdos.out'):
                cmd = "Yphon -tranI 2 -DebCut 0.5 " + " <superfij.out"
                print(cmd, " at ", voldir)
                output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)

            if len(self.Flat)==0:
                print ("Calling yphon to get f_vib, s_vib, cv_vib at ", phdir)
            with open("vdos.out", "r") as fp:
                f_vib, U_ph, s_vib, cv_vib, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT \
                    = ywpyphon.vibrational_contributions(self.T_vib, dos_input=fp, energyunit='eV')

            self.Flat.append(f_vib)
            self.Slat.append(s_vib)
            self.Clat.append(cv_vib)
            os.chdir( cwd )

        if len(self.Vlat)<=0:
            print("\nFATAL ERROR! cannot find required data from phonon collection for metadata tag:", self.tag,"\n")
            sys.exit()
        self.Slat = np.array(sort_x_by_y(self.Slat, self.Vlat))
        self.Clat = np.array(sort_x_by_y(self.Clat, self.Vlat))
        self.Flat = np.array(sort_x_by_y(self.Flat, self.Vlat))
        self.Vlat = np.array(sort_x_by_y(self.Vlat, self.Vlat))
        self.Dlat = np.full((len(self.Vlat)), 400.)
        self.volT = np.zeros(len(self.T_vib))
        self.GibT = np.zeros(len(self.T_vib))


    def check_vol(self):
        print ("\nChecking compatibility between qha/Yphon data and static calculation:\n")
        for i,v in enumerate (self.Vlat):
            print (v, self.energies[list(self.volumes).index(v)])
        if len(self.Vlat)!=len(self.volumes):
            print("\nWarning! The static/qha calculations are not inconsistent! Let me see if I can resolve it\n")
        _volumes = list(self.volumes)
        for i, vol in enumerate(_volumes):
            if vol not in self.Vlat:
                del _volumes[i]
                del self.energies[i]
                del self.dos_objs[i]
                print ("data in static calculation with volume=", vol, "is discarded")
        if len(_volumes)<len(self.volumes) : self.volumes = np.array(_volumes)
        if len(self.Vlat)==len(self.volumes) and len(self.Vlat)>=5:
            print("\nOK, I found all needed  data\n")
        else:
            print("xxxxxx", self.Vlat, self.volumes)
            print("\nFATAL ERROR! It appears that the calculations are not all done!\n")
            sys.exit()


    # get the energies, volumes and DOS objects by searching for the tag
    def find_static_calculations(self):
        static_calculations = self.vasp_db.collection.find({'$and':[ {'metadata.tag': self.tag}, {'adopted': True} ]})
        energies = []
        volumes = []
        dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
        structure = None  # single Structure for QHA calculation
        for calc in static_calculations:
            vol = calc['output']['structure']['lattice']['volume']
            #if vol_within (vol, volumes): 
            if vol in volumes: 
                print ("WARNING: skipped volume =", vol)
                continue
            volumes.append(vol)
            energies.append(calc['output']['energy'])
            dos_objs.append(self.vasp_db.get_dos(calc['task_id']))
            # get a Structure. We only need one for the masses and number of atoms in the unit cell.
            if structure is None:
                structure = Structure.from_dict(calc['output']['structure'])
                print(structure)
                print ("\n")
                self.formula_pretty = structure.composition.reduced_formula
                self.natoms = len(structure.sites)
                sa = SpacegroupAnalyzer(structure)
                self.phase = sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
    
        # sort everything in volume order
        # note that we are doing volume last because it is the thing we are sorting by!
    
        from dfttk.utils import sort_x_by_y
        self.energies = sort_x_by_y(energies, volumes)
        self.dos_objs = sort_x_by_y(dos_objs, volumes)
        #self.volumes = sort_x_by_y(volumes,volumes)
        self.volumes = np.array(list(map(float,sorted(volumes))))
        print ("found volumes from static calculations:", volumes)

        if self.phasename is None: self.phasename = self.formula_pretty+'_'+self.phase
        if not os.path.exists(self.phasename):
            os.mkdir(self.phasename)


    def get_static_calculations(self):
        t0 = min(self.T)
        t1 = max(self.T)
        td = (t1-t0)/(len(self.T)-1)
        if self.td < 0: td = self.td
        #print("xxxxx=", t0,t1,td)
        #theall = np.empty([len(prp_T), int((t1-t0)/td+1.5), len(self.volumes)])
        self.theall = np.empty([14, len(self.T), len(self.volumes)])
        for i,dos in enumerate(self.dos_objs):
            #print ("processing dos object at volume: ", self.volumes[i], " with nT =", len(self.T))
            prp_vol = runthelec(t0, t1, td, self.xdn, self.xup, self.dope, self.ndosmx, 
                self.gaussian, self.natom, dos=dos, _T=self.T, fout=sys.stdout, vol=self.volumes[i])
            """
            if 1==1:
                iFunc = trapz
            else:
                iFunc = simps
            F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=dos, fout=fvib, vol=volumes[i], IntegrationFunc=iFunc)
            Tuple of 14 float array containing
            thermal electron free energy, entropy, specific heat, M_el, seebeck_coefficients, 
            effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, temperature.
            Other quantities are for researching purpose
            """
            self.theall[:,:,i] = np.array( prp_vol ) # convert Tuple into array for the convenience of quasistatic interpolation


    def find_vibrational(self):
        if self.pyphon:
            self.toYphon()
            return

        if self.qhamode=="debye":
            self.qha_items = self.vasp_db.db['qha'].find({'metadata.tag': self.tag})
        elif self.qhamode=="phonon":
            self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata.tag': self.tag})
        else:
            try:
                self.qhamode='phonon'
                self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata.tag': self.tag})
            except:
                self.qhamode='debye'
                self.qha_items = self.vasp_db.db['qha'].find({'metadata.tag': self.tag})
    
        try:
            sys.stdout.write("\nTrying to get quasiharmonic mode : {}... \n".format(self.qhamode))
            self.T_vib = self.qha_items[0][self.qhamode]['temperatures'][::self.everyT]
        except:
            try:
                self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata.tag': self.tag})
                self.T_vib = self.qha_items[0][self.qhamode]['temperatures'][::self.everyT]
            except:
                self.toYphon()
                self.pyphon = True


    def get_qha(self):
        _Vlat = self.qha_items[0][self.qhamode]['volumes']
        _Slat = self.qha_items[0][self.qhamode]['entropies']
        _Clat = self.qha_items[0][self.qhamode]['heat_capacities']
        _Flat = self.qha_items[0][self.qhamode]['helmholtz_energies']
        _Dlat = self.qha_items[0]['debye']['debye_temperatures']
        Vlat = []
        Slat = []
        Clat = []
        Flat = []
        Dlat = []
        for i, vol in enumerate(_Vlat):
            if vol in Vlat: continue
            Vlat.append(vol)
            Slat.append(_Slat[i])
            Clat.append(_Clat[i])
            Flat.append(_Flat[i])
            Dlat.append(_Dlat[i])

        Slat = sort_x_by_y(Slat, Vlat)
        Clat = sort_x_by_y(Clat, Vlat)
        Flat = sort_x_by_y(Flat, Vlat)
        Dlat = sort_x_by_y(Dlat, Vlat)
        Vlat = sort_x_by_y(Vlat, Vlat)
        self.Slat = np.array(Slat)[:,::self.everyT]
        self.Clat = np.array(Clat)[:,::self.everyT]
        self.Flat = np.array(Flat)[:,::self.everyT]
        self.Dlat = np.array(Dlat)
        self.Vlat = np.array(Vlat)
        self.volT = self.qha_items[0][self.qhamode]['optimum_volumes'][::self.everyT]
        self.GibT = self.qha_items[0][self.qhamode]['gibbs_free_energy'][::self.everyT]
        if self.td < 0:
            s = []
            f = []
            c = []
            for i in range(len(self.volumes)):
                s.append(interp1d(self.T_vib, self.Slat[i,:])(self.T))
                f.append(interp1d(self.T_vib, self.Flat[i,:])(self.T))
                c.append(interp1d(self.T_vib, self.Clat[i,:])(self.T))
            self.Slat = np.array(s)
            self.Clat = np.array(c)
            self.Flat = np.array(f)
            self.volT = interp1d(self.T_vib, self.volT)(self.T)
            self.GibT = interp1d(self.T_vib, self.GibT)(self.T)
        self.Slat[np.isnan(self.Slat)] = 0

    
    def calc_thermal_expansion(self,i,kind='cubic'):
        if self.T[i]==0.0: beta=0.0
        if self.smooth:
            if self.eqmode==5:
                E0, Flat, Fel, Slat, Sel = BMsmooth(self.volumes, self.energies, self.Flat[:,i], 
                    self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:], BMvol5)
            else:
                #print (self.theall[0,i,:].shape, len(self.volumes), len(self.energies))
                E0, Flat, Fel, Slat, Sel = BMsmooth(self.volumes, self.energies, self.Flat[:,i], 
                    self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:], BMvol4)
        else:
            E0, Flat, Fel, Slat, Sel = self.energies, self.Flat[:,i], \
                self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:]

        #if True:
        try:
            if self.eqmode==4:
                blat, self.volT[i], self.GibT[i], P = BMDifB(self.volT[i], self.volumes, 
                    E0+Flat+Fel, BMvol4, N=7)
                if self.T[i]!=0.0: beta = (BMfitP(self.volT[i], self.volumes, E0+Flat+Fel+self.T[i]*(Slat+Sel), 
                        BMvol4) + P)/self.T[i]/blat
            elif self.eqmode==5:
                blat, self.volT[i], self.GibT[i], P = BMDifB(self.volT[i], self.volumes, 
                    E0+Flat+Fel, BMvol5, N=7)
                if self.T[i]!=0.0: beta = (BMfitP(self.volT[i], self.volumes, E0+Flat+Fel+self.T[i]*(Slat+Sel),
                        BMvol5) + P)/self.T[i]/blat
            else:
                blat, self.volT[i], newF = CenDifB(self.volT[i], self.volumes, 
                    E0+Flat+Fel, BMvol4, N=7,kind=kind)
                if newF!=None: self.GibT[i] = newF
                if self.T[i]!=0.0: beta = CenDif(self.volT[i], self.volumes, Slat+Sel, N=7,kind=kind)/blat
            return blat, beta
        except:
            return -1.0, 0.0

    
    def calc_thermodynamics(self):
        Faraday_constant = physical_constants["Faraday constant"][0]
        electron_volt = physical_constants["electron volt"][0]
        angstrom = 1e-30
        toJmol = Faraday_constant/self.natoms
        toGPa = electron_volt/angstrom*1.e-9
        thermofile = self.phasename+'/'+self.outf
        
        with open(thermofile, 'w') as fvib:
            fvib.write('#Found quasiharmonic mode : {}\n'.format(self.qhamode))
            if self.hasSCF:
                fvib.write('#T(K), volume, F(eV), S(J/K), H(J/K), a(-6/K), Cp(J/mol), Cv, Cpion, Bt(GPa), T_ph-D(K), T-D(K), F_el_atom, S_el_atom, C_el_atom, M_el, Seebeck_coefficients(μV/K), Lorenz_number(WΩK^{−2}), Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
            else:
                fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e Vol Gibbs_energy\n')
            for i in range(len(self.T)):
                if self.hasSCF:
                    if self.elmode==0: 
                        blat, beta = self.calc_thermal_expansion(i, kind='cubic')
                    else:
                        blat, beta = self.calc_thermal_expansion(i, kind='UnivariateSpline')
                    if blat < 0: 
                        print ("\nblat<0! Perhaps it has reached the upvolume limit at T =", self.T[i], "\n")
                        break
                    try:
                        slat = interp1d(self.volumes, self.Slat[:,i])(self.volT[i])
                        clat = interp1d(self.volumes, self.Clat[:,i])(self.volT[i])
                        flat = interp1d(self.volumes, self.Flat[:,i])(self.volT[i])
                        dlat = interp1d(self.volumes, self.Dlat)(self.volT[i])
                        cplat = clat+beta*beta*blat*self.volT[i]*self.T[i]
                    except:
                        print ("\nPerhaps it has reached the upvolume limit at T =", self.T[i], "\n")
                        break

                #print (self.theall.shape)
                prp_T = np.zeros((self.theall.shape[0]))
                for j in range(len(prp_T)):
                    if self.elmode==0: 
                        prp_T[j] = interp1d(self.volumes, self.theall[j,i,:], kind='cubic')(self.volT[i])
                    else:
                        prp_T[j] = UnivariateSpline(self.volumes, self.theall[j,i,:])(self.volT[i])
    
                # 0 K limit for the Lorenz number
                L = 2.443004551768e-08 #(1.380649e-23/1.60217662e-19*3.14159265359)**2/3
                if prp_T[5] > 1.e-16: L = prp_T[2]/prp_T[5]*k_B #avoiding divided by zero

                if self.hasSCF:
                    debyeT = get_debye_T_from_phonon_Cv(self.T[i], clat, dlat, self.natoms)
                    fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                    format(self.T[i], self.volT[i]/self.natoms, self.GibT[i]/self.natoms, (slat+prp_T[1])*toJmol,
                    (self.GibT[i]+self.T[i]*(slat+prp_T[1]))*toJmol, 
                    beta/3., (cplat+prp_T[2])*toJmol, (clat+prp_T[2])*toJmol, 
                    cplat*toJmol, blat*toGPa, debyeT, dlat, 
                    prp_T[0]/self.natoms, prp_T[1]*toJmol, prp_T[2]*toJmol, 
                    prp_T[3], prp_T[4], L, prp_T[5], prp_T[6], 
                    prp_T[7], prp_T[8], prp_T[10], prp_T[11], prp_T[12], prp_T[13]))
                else:
                    #(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, 
                    #Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i])
                    fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                    format(self.T[i], prp_T[0], prp_T[1], prp_T[2], prp_T[3], prp_T[4], L, 
                    prp_T[5], prp_T[6], prp_T[7], prp_T[8], prp_T[10], prp_T[11], prp_T[12], prp_T[13], 
                    self.volT[i], self.GibT[i]))
        return np.array(self.volumes)/self.natoms, np.array(self.energies)/self.natoms, thermofile


    def run_console(self):
        self.find_static_calculations()
        if not self.renew:
            pdis298 = self.phasename+'/figures/vdis298.15.png'
            if os.path.exists(pdis298): return None, None, None
      
        self.find_vibrational()

        T = np.array(self.T_vib)
        if self.pyphon:
            self.T = self.T_vib
        elif self.td < 0:
            self.T = T_remesh(min(self.T_vib), min(self.t1,max(self.T_vib)), self.td)
            #print ("xxxxx 2", len(self.T))
        else:
            self.T = T[T<=self.t1]
            #print ("xxxxx 3", len(self.T))

        if not self.pyphon: self.get_qha()
        self.hasSCF = True
        self.check_vol()
        if self.noel : self.theall = np.zeros([14, len(self.T), len(self.volumes)])
        else : self.get_static_calculations()
        return self.calc_thermodynamics()


if __name__ == '__main__':
    # initialize temperatures
    t0 = 0  # low temperature
    t1 = 1300  # high temperature
    td = 10  #
    # both double precision
    xdn = -100  # what is this XXX
    xup = 100  # what is this XXX
    # back to reals
    ndosmx = 10001  # 10001 #aka n_dos
    dope = 0.0
    natom = 1
    gaussian = 1000. #densed mesh near Fermi energy
    outf = "fvib_ele"

    # handling the command line option
    # TODO: use proper argparse module for this
    count = 1
    while (count < len(sys.argv)):
      if (sys.argv[count] == "-T1"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        t1 = float(sys.argv[count])
      elif (sys.argv[count] == "-dT"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        td = float(sys.argv[count])
      elif (sys.argv[count] == "-dope"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        dope = float(sys.argv[count])
        if abs(dope)<5.e-9:
          ndosmx = 100001
          gaussian = 10000.
      elif (sys.argv[count] == "-gauss"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        gaussian = float(sys.argv[count])
      elif (sys.argv[count] == "-ne"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        ndosmx = int(sys.argv[count])
      elif (sys.argv[count] == "-natom"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        natom = int(sys.argv[count])
      elif (sys.argv[count] == "-outf"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        outf = str(sys.argv[count])
      count = count + 1

    with open(outf, 'w') as fvib:
        F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom)
        fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
        for i in range(T.size):
            L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
            if Q_el[i] > 1.e-16: L = C_el_atom[i]/Q_el[i]*k_B
            fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i]))

