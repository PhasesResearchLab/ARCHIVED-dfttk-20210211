# -------- energy in eV, temperature in K
# assume every variable starting with a-h  and o-z are real numbers
# common block named comcon
from __future__ import division
import sys
import math
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d
from scipy.integrate import quadrature


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
    NELECTRONS = ados[idx] - e[idx]/(e[idx+1]-e[idx])*(ados[idx+1]-ados[idx])

    dF = 0.0
    #print ("NELECTRONS=", NELECTRONS)
    if dope != 0.0:
        NELECTRONS = NELECTRONS+dope
        idx = closest(ados,NELECTRONS)
        if idx == (NEDOS-1) or ados[idx] == ados[idx+1]:
            # we are dopidxng too much
            raise ValueError('Too much doping')
        dF = (NELECTRONS-ados[idx])/(ados[idx+1] - ados[idx])*(e[idx+1] - e[idx])+e[idx]
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

def gfind(mu_el, pe, pdos, NELECTRONS, Beta):
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
    return trapz(fn, pe[0:k])- NELECTRONS


# line 363
def caclf(pe, pdos, NELECTRONS, Beta, mu_ref=0.0, dF=0.0): #line 363
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
    mu_el = brentq(gfind, mu_ref-5.0, mu_ref+5.0, args=(pe, pdos, NELECTRONS, Beta), maxiter=10000)
    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k1 = len(tc)
    tf = 1.0/(np.exp(tc)+1.0)
    fn = pdos[0:k1]*pe[0:k1]*tf
    u = trapz(fn, pe[0:k1])

    k0 = closest(tc,-200)
    tf0 = tf[k0:]
    pdos = pdos[k0:k1]
    pe = pe[k0:k1]
    tf1 = 1.0 - tf0 + 1.e-60 # 1.e-60 is used to avoid log exception
    fn = pdos*(tf0*np.log(tf0)+tf1*np.log(tf1))
    s = trapz(fn, pe)

    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    fn2 = pdos*tf*(pe-mu_el)
    Q_el = trapz(fn, pe)
    Q_p = trapz(fn[pe<=dF], pe[pe<=dF])
    Q_e = trapz(fn[pe>dF], pe[pe>dF])
    Y_el = trapz(fn2, pe)

    fn = pdos*(pe-mu_el)*tf
    e_ = trapz(fn, pe)/Q_el
    fn = pdos[0:k1]*(pe[0:k1]-mu_el-e_)**2*tf
    cv = trapz(fn, pe[0:k1])
    fn = pdos[0:k1]*(pe[0:k1]-mu_el)**2*tf
    c_mu = trapz(fn, pe[0:k1])

#   hole/electron concentration by effective carrier
    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    f2 = interp1d(pe, fn, kind='linear')
    fmu = f2(mu_el)
    x = np.hstack([pe[pe<mu_el],mu_el])
    y = np.hstack([fn[pe<mu_el],fmu])
    W_p = trapz(y,x)
    x = np.hstack([mu_el, pe[pe>mu_el]])
    y = np.hstack([fmu, fn[pe>mu_el]])
    W_e = trapz(y,x)
    #W_e = trapz(fn[pe>mu_el], pe[pe>mu_el])
    #W_e = trapz(fn[pe>dF], pe[pe>dF])

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
    Y_p = trapz(y,x)

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
    Y_e = trapz(y,x)

    return mu_el, u, -s*k_B, cv*k_B*Beta*Beta, Q_el, Y_el, Q_p, Q_e, c_mu*k_B*Beta*Beta, W_p, W_e, Y_p, Y_e


def runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=sys.stdin, fout=sys.stdout, vol=None):
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
    T = np.arange(t0,t1+td,td) # temperature
    nT = T.size
    U_el = np.zeros(nT)
    S_el = np.zeros(nT)
    C_el = np.zeros(nT) # electronic specific heat
    C_mu = np.zeros(nT) # electronic specific heat at constant chemical potential
    M_el = np.zeros(nT) # electronic chemical potential
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
        M_el[i], U_el[i], S_el[i], C_el[i], Q_el[i],Y_el[i], Q_p[i],Q_e[i],  C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i] = caclf(e, dos, NELECTRONS, Beta, M_el[i-1], dF=-dF)
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


def CenDif(v, vol, F, N=7,kind='cubic'):
    dV = 0.001*(max(vol) - min(vol))
    if kind=='cubic':
        f2 = interp1d(vol, F, kind='cubic')
    else:
        from scipy.interpolate import UnivariateSpline
        f2 = UnivariateSpline(vol, F)
    result = 0.0 
    n = N//2
    for i in range(0, n):
        #print (v, f2(v+dV*(i+1)) , f2(v-dV*(i+1)))
        result += (f2(v+dV*(i+1)) - f2(v-dV*(i+1)))/(2.0*(i+1)*dV)
    #sys.exit()
    return result/n


def CenDifB(v, vol, F, N=7,kind='cubic'):
    dV = 0.001*(max(vol) - min(vol))
    result = 0.0 
    n = N//2
    for i in range(0, n):
        result += (CenDif(v+dV*(i+1), vol, F, kind=kind) - CenDif(v-dV*(i+1), vol, F, kind=kind))/(2.0*(i+1)*dV)
    return result*v/(n)


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

def get_debye_T_from_phonon_Cv(temperature, clat, dlat, natoms, _td=50):
    t0 = dlat
    d0 = debye_heat_capacity(temperature, t0, natoms) - clat
    #print("lat=", temperature, 96484*debye_heat_capacity(temperature, t0, natoms)/natoms, t0, 96484*clat/natoms)
    if d0 > 0.0: td = _td
    elif d0 <0.0: td = -_td 
    else: return t0
    for i in range(999):
        if t0 < 0.1 : return 0
        t1 = t0 + td
        t1 = max(t1,0.01)
        d1 = debye_heat_capacity(temperature, t1, natoms) - clat
        if d1*d0 < 0.0: break
        t0 = t1
        d0 = d1
        
    for i in range(20):
        if abs(t1-t0) < 0.0001: break
        t2 = 0.5*(t0 + t1)
        d2 = debye_heat_capacity(temperature, t2, natoms) - clat
        if d2*d0 < 0.0:
            t1 = t2
            d1 = d2
        elif d2*d0 > 0.0:
            t0 = t2
            d0 = d2
        else:
            return t2

    return 0.5*(t0 + t1)

def thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, db_file, everyT=1, metadata=None, qhamode=None, plot=False):
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
    metadata : str
        metadata tag to access the to be calculated the compound
    qhamode : str
        Mode for the quasiharmonic calculations, according to it the T dependence volume to be extracted

    Output to (printed to outf)
    ---------------------------
    The properties in the order of temperature, thermal electron free energy, entropy, 
    specific heat, seebeck_coefficients, Lorenz number,  
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, ..., 
    and the equilibrium volume extracted from MongoDB in the last column
    """

    from atomate.vasp.database import VaspCalcDb
    from pymatgen import Structure
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)

    # get the energies, volumes and DOS objects by searching for the tag
    static_calculations = vasp_db.collection.find({'$and':[ {'metadata.tag': metadata}, {'adopted': True} ]})

    energies = []
    volumes = []
    dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
    structure = None  # single Structure for QHA calculation
    for calc in static_calculations:
        vol = calc['output']['structure']['lattice']['volume']
        if vol in volumes: continue
        volumes.append(vol)
        energies.append(calc['output']['energy'])
        dos_objs.append(vasp_db.get_dos(calc['task_id']))
        # get a Structure. We only need one for the masses and number of atoms in the unit cell.
        if structure is None:
            structure = Structure.from_dict(calc['output']['structure'])
            print(structure)

    # sort everything in volume order
    # note that we are doing volume last because it is the thing we are sorting by!

    from dfttk.utils import sort_x_by_y
    energies = sort_x_by_y(energies, volumes)
    dos_objs = sort_x_by_y(dos_objs, volumes)
    volumes = np.array(list(map(float,sorted(volumes))))
    if qhamode=="debye":
        qha_items = vasp_db.db['qha'].find({'metadata.tag': metadata})
    elif qhamode=="phonon":
        qha_items = vasp_db.db['qha_phonon'].find({'metadata.tag': metadata})
    else:
        try:
            qhamode='phonon'
            qha_items = vasp_db.db['qha_phonon'].find({'metadata.tag': metadata})
        except:
            qhamode='debye'
            qha_items = vasp_db.db['qha'].find({'metadata.tag': metadata})

    try:
        T = qha_items[0][qhamode]['temperatures'][::everyT]
    except:
        try:
            print("\nno data entry found for quasiharmonic mode :", qhamode, "\n")
            print("trying to get from phonon.\n")
            qha_items = vasp_db.db['qha_phonon'].find({'metadata.tag': metadata})
            T = qha_items[0][qhamode]['temperatures'][::everyT]
        except:
            print("\nERROR STOP! no data entry found for quasiharmonic mode :", qhamode, "\n")
            sys.exit()
    T = np.array(T)
    T= T[T<=t1]

    natoms = qha_items[0][qhamode]['natoms']
    eq_energy = qha_items[0]['eos']['eq_energy']
    import scipy
    Faraday_constant = scipy.constants.physical_constants["Faraday constant"][0]
    electron_volt = scipy.constants.physical_constants["electron volt"][0]
    angstrom = 1e-30
    toJmol = Faraday_constant/natoms
    toGPa = electron_volt/angstrom*1.e-9
    hasSCF = False
    #try:
    if True:
        Slat = qha_items[0][qhamode]['entropies']
        Clat = qha_items[0][qhamode]['heat_capacities']
        Flat = qha_items[0][qhamode]['helmholtz_energies']
        Vlat = qha_items[0][qhamode]['volumes']
        volumes = np.array(list(map(float,sorted(Vlat))))
        Slat = sort_x_by_y(Slat, volumes)
        Clat = sort_x_by_y(Clat, volumes)
        Flat = sort_x_by_y(Flat, volumes)
        Slat = np.array(Slat)[:,::everyT]
        Clat = np.array(Clat)[:,::everyT]
        Flat = np.array(Flat)[:,::everyT]
        Dlat = qha_items[0]['debye']['debye_temperatures']
        Dlat = sort_x_by_y(Dlat, volumes)
        Dlat = np.array(Dlat)
        hasSCF = True
    #except:
    #    pass

    print("\nFound quasiharmonic mode :", qhamode, "\n")
    volT = qha_items[0][qhamode]['optimum_volumes'][::everyT]
    GibT = qha_items[0][qhamode]['gibbs_free_energy'][::everyT]
    t0 = min(T)
    t1 = max(T)
    td = (t1-t0)/(len(T)-1)

    prp_T = np.empty([14])
    theall = np.empty([len(prp_T), int((t1-t0)/td+1.5), len(volumes)])
    with open(outf, 'w') as fvib:
        fvib.write('#Found quasiharmonic mode : {}\n'.format(qhamode))
        if hasSCF:
            fvib.write('#T(K), volome, F(eV), S(J/K), H(J/K), a(-6/K), Cp(J/mol), Cv, Cpion, T_ph-D(K), T-D(K),            Cv  Bt(GPa), F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
        else:
            fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e Vol Gibbs_energy\n')
        for i,dos in enumerate(dos_objs):
            print ("processing dos object at volume: ", volumes[i], " with nT =", int((t1-t0)/td+1.5))
            prp_vol = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=dos, fout=fvib, vol=volumes[i])
            """
            F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=dos, fout=fvib, vol=volumes[i])
            Tuple of 14 float array containing
            thermal electron free energy, entropy, specific heat, M_el, seebeck_coefficients, 
            effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, temperature.
            Other quantities are for researching purpose
            """
            if len(prp_T) != len(prp_vol):
                print("ERROR! Length of prp_T of ", len(prp_T)," is not equal to that of ", len(prp_vol), " for prp_vol")
                sys.exit()
            theall[:,:,i] = np.array( prp_vol ) # convert Tuple into array for the convenience of quasistatic interpolation
        #f2 = interp1d(T, volT, kind='cubic')
        #TT = np.linspace(min(T),max(T),10000)
        #VV = f2(TT)
        #DD = np.zeros(VV.shape,np.float)
        #DD[0:-1] = np.diff(VV) / np.diff(TT)
        #DD[-1] = (VV[-1] - VV[-2])/(TT[-1] - TT[-2])
        #DD /= VV
        for i in range(len(T)):
            for j in range(len(prp_T)):
                f2 = interp1d(volumes, theall[j,i,:], kind='cubic')
                prp_T[j] = f2(volT[i])
            if hasSCF:
                f2 = interp1d(volumes, Slat[:,i])
                slat = f2(volT[i])
                f2 = interp1d(volumes, Clat[:,i])
                clat = f2(volT[i])
                f2 = interp1d(volumes, Flat[:,i])
                flat = f2(volT[i])
                f2 = interp1d(volumes, Dlat)
                dlat = f2(volT[i])
                #blat = CenDifB(volT[i], volumes, energies+Flat[:,i]+theall[0,i,:],N=7,kind='cubic')
                #beta = CenDif(volT[i], volumes, Slat[:,i]+theall[1,i,:],N=7,kind='cubic')/blat
                blat = CenDifB(volT[i], volumes, energies+Flat[:,i]+theall[0,i,:],N=7,kind='UnivariateSpline')
                beta = CenDif(volT[i], volumes, Slat[:,i]+theall[1,i,:],N=7,kind='UnivariateSpline')/blat
                cplat = clat+beta*beta*volT[i]*T[i]
            #f2 = interp1d(TT, DD, kind='cubic')
            #beta = f2(T[i])
            L = 2.443004551768e-08 #(1.380649e-23/1.60217662e-19*3.14159265359)**2/3 # 0 K limit for the Lorenz number
            if prp_T[5] > 1.e-16: L = prp_T[2]/prp_T[5]*k_B #avoiding divided by zero
            if hasSCF:
                debyeT = get_debye_T_from_phonon_Cv(T[i], clat, dlat, natoms)
                fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                format(T[i], volT[i]/natoms, GibT[i]/natoms, (slat+prp_T[1])*toJmol,
                (GibT[i]+T[i]*(slat+prp_T[1]))*toJmol, 
                beta/3., (cplat+prp_T[2])*toJmol, (clat+prp_T[2])*toJmol, 
                cplat*toJmol, blat*toGPa, debyeT, dlat, 
                prp_T[0]/natoms, prp_T[1]*toJmol, prp_T[2]*toJmol, 
                prp_T[3], prp_T[4], L, prp_T[5], prp_T[6], 
                prp_T[7], prp_T[8], prp_T[10], prp_T[11], prp_T[12], prp_T[13]))
                #(GibT[i]+T[i]*(slat+prp_T[1])-eq_energy)*toJmol, 
            else:
                fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], prp_T[0], prp_T[1], prp_T[2], prp_T[3], prp_T[4], L, prp_T[5], prp_T[6], prp_T[7], prp_T[8], prp_T[10], prp_T[11], prp_T[12], prp_T[13], volT[i], GibT[i]))
            #fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i]))
   
    return np.array(volumes)/natom, np.array(energies)/natoms


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
    metadata = None

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

