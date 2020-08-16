#!/global/project/projectdirs/m891/yiwang62/anaconda3/bin/python
#!/usr/bin/python -x

import sys
from datetime import datetime

import os, fnmatch
import copy
import time
import datetime
import numpy as np
from scipy.optimize import linprog
from scipy.interpolate import interp1d
from scipy import interpolate
from numpy.linalg import solve
from fractions import Fraction
from difflib import SequenceMatcher

#from PIL import Image
from scipy import misc
from scipy import ndimage as ndi
import math
import glob
from scipy.optimize import curve_fit
from scipy.constants import physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz, simps

import re
import json
import subprocess
from shutil import copyfile

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from difflib import SequenceMatcher
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from elements import elements

MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0}

periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""

def SGTE(T,a):
  fval = a[0]+a[1]*T
  if len(a) > 2:
    fval += a[2]*T*np.log(T)
  if len(a) > 3:
    fval += a[3]*T*T
  if len(a) > 4:
    fval += a[4]*T*T*T
  if len(a) > 5:
    fval += a[5]/T
  return(fval)

def SGTE2(T, a, b):
  return (SGTE(T, [a,b]))

def SGTE3(T, a, b, c):
  return (SGTE(T, [a,b,c]))

def SGTE4(T, a, b, c, d):
  return (SGTE(T, [a,b,c,d]))

def SGTE5(T, a, b, c, d, e):
  return (SGTE(T, [a,b,c,d,e]))

def SGTE6(T, a, b, c, d, e, f):
  return (SGTE(T, [a,b,c,d,e,f]))


def SGTEC1(T,a):
  return C_SGTE(T,[a])

def SGTEC2(T,a,b):
  return C_SGTE(T,[a,b])

def SGTEC3(T,a,b,c):
  return C_SGTE(T,[a,b,c])

def SGTEC4(T,a,b,c,d):
  return C_SGTE(T,[a,b,c,d])

def C_SGTE(T,a):
  fval = 0
  if len(a) > 0:
    fval += a[0]
  if len(a) > 1:
    fval += a[1]*T
  if len(a) > 2:
    fval += a[2]*T*T
  if len(a) > 3:
    fval += a[3]/T/T
  return(fval)


def SGTES(T,f):
  s = 0.0
  if len(f)>0:
    s += f[0]
  if len(f)>1:
    s += f[1]*np.log(T)
  if len(f)>2:
    s += f[2]*T
  if len(f)>3:
    s += f[3]*T*T
  if len(f)>4:
    s += f[4]/T/T
  return s

def SGTEH(T,f):
  h = 0.0
  if len(f)>0:
    h += f[0]
  if len(f)>1:
    h += f[1]*T
  if len(f)>2:
    h += f[2]*T*T
  if len(f)>3:
    h += f[3]*T*T*T
  if len(f)>4:
    h += f[4]/T
  return h

def SGTEC(T,f):
  s = 0.0
  if len(f)>0:
    s += f[0]
  if len(f)>1:
    s += f[1]*T
  if len(f)>2:
    s += f[2]*T*T
  if len(f)>3:
    s += f[3]/T/T
  return s


def CSGTEfit(f, x, y):
  popt,pcov = curve_fit(f, x, y)
  z = C_SGTE(x,popt)
  ferror=math.sqrt(((z-y)**2).sum()/len(z))
  return(popt,ferror)


def fitStoichiometricCp(x,y, thr=0.001):
  f,ferror = CSGTEfit(SGTEC2, x, y)
  #if ferror > thr:
  #  f,ferror = CSGTEfit(SGTEC2, x, y)
  if ferror > thr:
    f,ferror = CSGTEfit(SGTEC3, x, y)
  if ferror > thr:
    f,ferror = CSGTEfit(SGTEC4, x, y)
  return f,ferror


def H_SGTE(T,c):
  h = 0.
  if len(c)>0:
    h += c[0]*T
  if len(c)>1:
    h += c[1]/2*T*T
  if len(c)>2:
    h += c[2]/3*T*T*T
  if len(c)>3:
    h += -c[3]/T
  return h

def fitStoichiometricH(x,y,c):
  zz = H_SGTE(x,c)
  h = (y - zz).sum()/len(y)
  ferror=math.sqrt(((h+zz-y)**2).sum()/len(zz))
  h = [h]
  if len(c)>0:
    h.append(c[0])
  if len(c)>1:
    h.append(c[1]/2)
  if len(c)>2:
    h.append(c[2]/3)
  if len(c)>3:
    h.append(-c[3])
  return h,ferror


def S_SGTE(T,c):
  s = 0.
  if len(c)>0:
    s += c[0]+c[0]*np.log(T)
  if len(c)>1:
    s += c[1]*T
  if len(c)>2:
    s += c[2]/2*T*T
  if len(c)>3:
    s += -c[3]/2/T/T
  return s

def fitStoichiometricS(x,y,c):
  zz = S_SGTE(x,c)
  b = (y - zz).sum()/len(y)
  ferror=math.sqrt(((b+zz-y)**2).sum()/len(zz))
  s = []
  if len(c)>0:
    s.append(b+c[0])
    s.append(c[0])
  if len(c)>1:
    s.append(c[1])
  if len(c)>2:
    s.append(c[2]/2)
  if len(c)>3:
    s.append(-c[3]/2)
  return s,ferror

def fitStoichiometric(x,y, thr=1.0):
  f,ferror = SGTEfit(SGTE2, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE3, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE4, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE5, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE6, x, y)
  return f,ferror

def SGTEfit(f, x, y):
  popt,pcov = curve_fit(f, x, y)
  z = SGTE(x,popt)
  ferror=math.sqrt(((z-y)**2).sum()/len(z))
  return(popt,ferror)

def outexpressionG(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*log(T)'.format(f)
      elif i==3:
        out += ' {:+g}*T*T'.format(f)
      elif i==4:
        out += ' {:+g}*T*T*T'.format(f)
      elif i==5:
        out += ' {:+g}/T'.format(f)
    return out 

def outexpressionS(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*log(T)'.format(f)
      elif i==2:
        out += ' {:+g}*T'.format(f)
      elif i==3:
        out += ' {:+g}*T*T'.format(f)
      elif i==4:
        out += ' {:+g}/T/T'.format(f)
    return out

def outexpressionH(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*T'.format(f)
      elif i==3:
        out += ' {:+g}*T*T*T'.format(f)
      elif i==4:
        out += ' {:+g}/T'.format(f)
    return out

def outexpressionCp(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*T'.format(f)
      elif i==3:
        out += ' {:+g}/T/T'.format(f)
    return out


def proStoichiometricG():
    #try:
      x = zthermo.get("temperature (K)")
      y = zthermo.get("Gibbs energy (eV/atom)")
      H298 = threcord.get("H298.15 (J/mol-atom)")
      x = np.array(list(map(float, x)))
      y = np.array(list(map(float, y)))*eVtoJ - H298
      i0 = 0
      for i,T in enumerate(x):
        if T < T0:
          i0 = i
      ifit0 = i0-15
      ifit0 = max(ifit0,0)
      
      f,ferror = fitStoichiometric(x[ifit0:],y[ifit0:])
      gout = 'G(T) =' + outexpressionG(f)
      print("fitting uncertainty=", ferror)
      #print(gout)
      s = []
      h = []
      c = []
      if len(f) >0:
        h.append(f[0])
      if len(f) >1:
        s.append(-f[1])
      if len(f) >2:
        s = []
        s.append(-f[1]-f[2])
        s.append(-f[2])
        h.append(-f[2])
        c.append(-f[2])
      if len(f) >3:
        s.append(-2.0*f[3])
        h.append(-f[3])
        c.append(-2.0*f[3])
      if len(f) >4:
        s.append(-3.0*f[4])
        h.append(-2.0*f[4])
        c.append(-6.0*f[4])
      if len(f) >5:
        s.append(f[5])
        h.append(2.0*f[5])
        c.append(-2.0*f[5])
      sout = 'S(T) =' + outexpressionS(s)
      hout = 'H(T) =' + outexpressionH(h)
      cout = 'Cp(T) =' + outexpressionCp(c)
      """
      print (sout)
      print (hout)
      print (cout)
      """
      uncertanty = {}
      SGTErec.update({"G-H298.15 (J/mol-atom)":gout})
      SGTErec.update({"H-H298.15 (J/mol-atom)":hout})
      SGTErec.update({"S (J/mol-atom/K)":sout})
      SGTErec.update({"Cp (J/mol-atom/K)":cout})
      return(f,h,s,c,x[i0:])

    #except:
      pass


def proStoichiometricCp():
    #try:
      uncertanty = {}
      x = zthermo.get("temperature (K)")
      y = zthermo.get("Cp (J/mol-atom/K)")
      H298 = threcord.get("H298.15 (J/mol-atom)")
      x = np.array(list(map(float, x)))
      y = np.array(list(map(float, y)))
      i0 = 0
      for i,T in enumerate(x):
        if T < T0:
          i0 = i
      ifit0 = i0
      ifit0 = max(ifit0,0)
      #print("xxxxxxxx=",x[ifit0:],T0) 
      c,cerror = fitStoichiometricCp(x[ifit0:],y[ifit0:])

      y = zthermo.get("enthalpy (J/mol-atom)")
      y = np.array(list(map(float, y))) - H298
      h,herror = fitStoichiometricH(x[ifit0:],y[ifit0:],c)

      y = zthermo.get("entropy (J/mol-atom/K)")
      y = np.array(list(map(float, y)))
      s,serror = fitStoichiometricS(x[ifit0:],y[ifit0:],c)

      f = [h[0]]
      if len(s) >0:
        f.append(-s[0]+c[0])
        f.append(-c[0])
      if len(c) >1:
        f.append(-c[1]/2)
      if len(c) >2:
        f.append(-c[2]/6)
      if len(c) >3:
        f.append(-c[3]/2)
      gout = 'G(T) =' + outexpressionG(f)
      #print (gout)

      SGTErec.update({"T":[x[ifit0], x[-1]]})
      SGTErec.update({"Cp (J/mol-atom/K)":[outexpressionCp(c),{"error":round(cerror,2)}]})
      SGTErec.update({"H-H298.15 (J/mol-atom)":[outexpressionH(h),{"error":round(herror,2)}]})
      SGTErec.update({"S (J/mol-atom/K)":[outexpressionS(s),{"error":round(serror,2)}]})
      SGTErec.update({"G-H298.15 (J/mol-atom)":[outexpressionG(f),{"error":round(herror,2)}]})
      return(f,h,s,c,x[i0:])

    #except:
      pass

def thermoplot(folder,thermodynamicproperty,x,y,yzero=None,fitted=None,xT=None,xlabel="T (K)", xlim=None, elonly=None, expt=None, CoT=False):
    global mindex
    mindex = 0
    cwd = os.getcwd()
    os.chdir( folder )
    plt.rc('font', size=24)
    fig,ax=plt.subplots()
    fig.set_size_inches(12,9)
    ax.yaxis.set_ticks_position('both')
    if xlim!=None: ax.set_xlim([0,xlim])
    else: ax.set_xlim([0,np.array(list(map(float,x))).max()])
    fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+".png"
    if thermodynamicproperty.split('(')[0].strip()=="Debye temperature":
        if float(x[0])==0.0: y[0]=float("nan")

    if thermodynamicproperty=="0 K total energies (eV/atom)":
        #plt.xlabel("atomic volume (Angstrom^3)")
        xlabel = "atomic volume ($Angstrom^3$)"
        plt.ylabel("0 K total energies (eV/atom)")
        ax.set_xlim([min(x)*0.95,max(x)*1.05])
        ax.plot(x, y, fillstyle='none', marker='o', markersize=12, color='k', linestyle='None', label="DFT")
        xnew = np.linspace(min(x)*0.95,max(x)*1.05, 300)  
        #from scipy.interpolate import UnivariateSpline
        from dfttk.pythelec import BMvol4, BMvol
        f2, pcov = curve_fit(BMvol4, x, y)
        #f2 = UnivariateSpline(x,y)
        #f2 = interp1d(x, y)
        ynew = BMvol(xnew, f2)
        ax.plot(xnew,ynew,'-',linewidth=1,color='b', label="BMvol4")
    elif thermodynamicproperty!="heat capacities (J/mol-atom/K)":
      if yzero != None:
        y0 = np.nanmin(np.array(list(map(float,y))))
        y1 = np.nanmax(np.array(list(map(float,y))))
        ylow = 0.0
        yhigh = y1*1.05
        if y0 < 0.0:
          ylow = y0
        ax.set_ylim([ylow,yhigh])
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,4))
        #if yhigh > 1.e-2 : ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        #else : ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
      ax.plot(x,y,'-',linewidth=2,color='b', label=thermodynamicproperty)
      if fitted!=None:
        ax.plot(xT[::5],fitted[::5],'--',fillstyle='none', marker='o', markersize=12, linewidth=2,color='k', label="fitted")
      if xlim!=None: 
          t0 = np.array(list(map(float,x)))
          t1 = np.array(list(map(float,y)))
          tmp = np.array(list(map(float,t1[t0<=xlim])))
          ax.set_xlim([0.0,xlim])
          #print([0.95*np.nanmin(tmp),1.05*np.nanmax(tmp)])
          ax.set_ylim([0.98*np.nanmin(tmp),1.02*np.nanmax(tmp)])
          fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(xlim)+".png"
    else:
      if fitted!=None:
        y0 = []
        y1 = []
        y2 = []
        for x0,x1,x2 in y:
          y0.append(x0)
          y1.append(x1)
          y2.append(x2)
        ax.set_ylim([0.0,np.array(list(map(float,y0))).max()*1.05])
        ax.plot(x,y0,'-',linewidth=2,color='b', label="$C_p$")
        ax.plot(xT[::5],fitted[::5],'--',fillstyle='none', marker='o', markersize=12, linewidth=2,color='k', label="fitted")
        ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_v$")
        ax.plot(x,y2,':',linewidth=2,color='g', label="$C_{v,ion}$")
        y2 = np.array(list(map(float,y1))) - np.array(list(map(float,y2)))
        ax.plot(x,y2,'-.',linewidth=2,color='r', label="$C_{el}$")
      else:
        y0 = []
        y1 = []
        for x0,x1 in y:
          y0.append(x0)
          y1.append(x1)
        if xlim!=None: 
            ylim = 0.0
            if CoT:
                xx = np.array(x)
                yy = np.array(y0)
                yy = yy[xx!=0.0]
                xx = xx[xx!=0.0]
                yy = yy/xx
                #xnew = np.linspace(0, math.sqrt(xlim), 100)
                xnew = np.linspace(min(x), math.sqrt(xlim), 10)
                xnew = xnew*xnew
                xx = xx*xx
                from scipy.interpolate import UnivariateSpline
                f2 = UnivariateSpline(xx, yy, s=0)
                ynew = f2(xnew)
                ax.set_xlim([0.0,xlim])
                ax.set_ylim([0.0,max(ynew)*1.5])
                xnew = xx
                ynew = yy
                fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(xlim)+"_T2.png"
            else:
                for i,v in enumerate(x):
                    if v <= xlim: ylim=max(ylim,y0[i])
                ax.set_xlim([0.0,xlim])
                ax.set_ylim([0.0,ylim*1.05])
                fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(xlim)+".png"
        elif elonly!=None:
            ylim = 0.0
            for i,v in enumerate(x):
                tmp = 0.0
                if v!=0.0: tmp = (y0[i]-y1[i])/v
                if v <= elonly: 
                    if CoT:
                        if v>0.0: ylim=max(ylim,tmp)
                        else: break
                    else:
                        if v>0.0: ylim=max(ylim,(y0[i]-y1[i]))
            ax.set_xlim([0.0,elonly])
            if ylim==0.0: ylim=1.e-4
            ax.set_ylim([0.0,ylim*1.5])
            tmp = ""
            if CoT: tmp = "_T"
            fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(elonly)+'_el'+tmp+".png"
        else:
            fname = thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+"2.png"
            ax.set_ylim([0.0,np.array(list(map(float,y0))).max()*1.05])
        if elonly==None:
            if CoT:
                ax.plot(xnew,ynew,'-',linewidth=2,color='b', label="$C_{p,lat+el}/T$")
                plot_expt(expt, 'heat capacity', ax, CoT=CoT, xlim=xlim)
            else:
                ax.plot(x,y0,'-',linewidth=2,color='b', label="$C_{p,lat+el}$")
                ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_{p,lat}$")
                plot_expt(expt, 'heat capacity', ax)
        y2 = np.array(list(map(float,y0))) - np.array(list(map(float,y1)))
        if CoT:
            if xlim==None:
                for i,v in enumerate(x):
                    if v>0.0: y2[i] /= v
                    else: y2[i] = float("nan")
                ax.plot(x,y2,'-.',linewidth=2,color='r', label="$C_{el}/T$")
                plot_expt(expt, 'electronic heat capacity', ax, CoT=CoT)
        else:
            ax.plot(x,y2,'-.',linewidth=2,color='r', label="$C_{el}$")
            plot_expt(expt, 'electronic heat capacity', ax, CoT=CoT)

    if CoT and xlim!=None:
        plt.xlabel("$T^2 (K^2)$")
    else:
        plt.xlabel(xlabel)
    if CoT:
        plt.ylabel("$C/T$ (J/mol-atom/K/K)")
    else:
        plt.ylabel(thermodynamicproperty)
    ax.legend(loc=0, prop={'size': 24})
    figures.update({thermodynamicproperty:folder.split('/')[-1]+'/'+fname})
    fig.savefig(fname,bbox_inches='tight')
    plt.close(fig)
    os.chdir( cwd )

def plot_expt (expt, prp, ax, CoT=False, xlim=None):
    global mindex
    if expt!=None:
        for k1 in expt:
            val1 = expt[k1]
            for k2 in val1:
                if k2=='property':
                    if val1[k2]==prp:
                        xval = np.array(val1['T'])
                        yval = np.array(val1['val'])
                        Author = val1['Author']
                        Unit = val1['Unit']
                        natom = val1['natom']
                        yval /= natom
                        if Unit=='mJ/K' : yval /= 1000.
                        if CoT:
                            if xlim!=None:
                                ax.plot(xval*xval,yval/xval, marker=markers[mindex%len(markers)], markersize=8, linestyle='None', label=Author.split(',')[0])
                            else:
                                ax.plot(xval,yval/xval, marker=markers[mindex%len(markers)], markersize=8, linestyle='None', label=Author.split(',')[0])
                        else:
                            ax.plot(xval,yval, marker=markers[mindex%len(markers)], markersize=8, linestyle='None', label=Author.split(',')[0])
                        mindex += 1

def myjsonout(data,fp,indent="",comma=""):
	#print (data)
	mj = ''
	if (isinstance(data,dict)):
		fp.write('{}\n'.format('{'))
			#sys.stdout.write('\n{}{}\n'.format(indent, '{'))
		nkey = 0
		for key in sorted(set(data.keys())):
			nkey += 1
			if nkey!=len(data):
				comma1 = ","
			else:
				comma1 = ""
			val = data[key]
			jval = json.dumps(val)
			jkey = json.dumps(key)
			#print (val)
			if (isinstance(val,dict)):
				fp.write('{}{}: '.format(indent+"    ",jkey))
				myjsonout(val,fp,indent+"    ",comma1)
			elif (isinstance(val,tuple)):
				#print (val)
				out = list(val)
				#print(out)
				fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, out, comma1))
			elif (isinstance(val,str)):
				if (indent == ""):
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
			else:
				if (indent==""):
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))

				#print(val)
				"""
				if (nkey!=len(data)):
					sys.stdout.write('{}{}: {},\n'.format(indent+"    ", key, val))
				else:
					sys.stdout.write('{}{}: {}\n'.format(indent+"    ", key, val))
				"""
		if comma==',':
			fp.write('{}{}{}\n\n'.format(indent,'}', comma))
		else:
			fp.write('{}{}{}\n'.format(indent, '}', comma))

def similar(pp,pall):
  known = ["L12", "delta", "D022", "Gamma"]
  ii = -1
  for o in known:
    if pp.find(o)>-1:
      pname = o
      ii = 0
      break
  if ii == -1:
    return "unknown"

  s = 0.0
  for i,p in enumerate(pall):
    snew = SequenceMatcher ( None, pname, p ).ratio()
    if snew > s:
      ii = i
      s = snew
  print (pp, "= ", pall[ii], " by ", s)
  if s > 0.5:
    return pall[ii]
  else:
    return "unknown"


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def formula2composition(formula):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'')
  newc = ""
  """Follow the convention, elemental symbol must start from capital letter"""
  for c in formula:
    if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      newc = newc + '|'
    newc = newc + c
  els = newc.split('|')
  els = [k for k in els if k != '']

  """now get the composition for each element"""
  ele = []
  com = []
  for el in els:
    newel = ""
    newcc = ""
    for c in el:
      if c.isalpha():
        newel = newel + c
      else:
        newcc = newcc + c

    if (newel not in periodictable):
      print('"',newel,'" is not an element! your formula is wrong!')
      sys.exit(1)
    ele.append(newel)

    if (len(newcc)!=0):
      if (isfloat(newcc)):
        com.append(int(newcc))
      else:
        print('"',newcc,'" is not an int number! your formula is wrong!')
        sys.exit(1)
    else:
      com.append(1.0)
  com = np.array(list(map(int,com)))

  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=int)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]
  return elist,clist

def formula2elist(formula):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'')
  newc = ""
  """Follow the convention, elemental symbol must start from capital letter"""
  for c in formula:
    if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      newc = newc + '|'
    newc = newc + c
  els = newc.split('|')
  els = [k for k in els if k != '']

  """now get the composition for each element"""
  ele = []
  com = []
  for el in els:
    newel = ""
    newcc = ""
    for c in el:
      if c.isalpha():
        newel = newel + c
      else:
        newcc = newcc + c

    if (newel not in periodictable):
      print('"',newel,'" is not an element! your formula is wrong!')
      sys.exit(1)
    ele.append(newel)

    if (len(newcc)!=0):
      if (isfloat(newcc)):
        com.append(float(newcc))
      else:
        print('"',newcc,'" is not a float number! your formula is wrong!')
        sys.exit(1)
    else:
      com.append(1.0)
  com = np.array(list(map(float,com)))
  com = com/sum(com)
  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=float)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]

  return elist

def prety_formulaO(longphasename):
  puc = longphasename.split('|')[-1]
  _els,_nat=formula2composition(puc)

def prety_formula(_els,_nat):
  els = sorted(set(_els))
  nat = np.zeros(len(els),dtype=int)
  for i,el in enumerate(_els):
    ix = els.index(el)
    nat[ix] += _nat[i]

  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  form = ""
  for j,el in enumerate(els):
    ix = nat[j]//i
    form = form+el
    if ix!=1:
      form = form+str(ix)
  return form


def Genergy(thermofile,dir0):
  tmelt = 9999.
  ele = threcord.get("Elements")
  if ele!=None:
    if len(ele)==1:
      tmelt = ELEMENTS[ele[0]].tmelt

  folder = dir0+'/'+"figures"
  if not os.path.exists(folder):
    os.mkdir(folder)

  vdos_e_Cij = thermofile.replace('/vdos_e','/vdos_e_Cij')
  #print(" I am he",  vdos_e_Cij)
  if not os.path.exists(vdos_e_Cij) : 
    vdos_e_Cij = thermofile.replace('/vdos_e','/vdos_Cij')
  if os.path.exists(vdos_e_Cij) :
    vdos_e_Cij = np.loadtxt(vdos_e_Cij, comments="#", dtype=np.float)
    ij = 0
    for i in range(1,7):
      for j in range(i,7):
        ij = ij + 2
        if abs(vdos_e_Cij[:,ij]).max() > 1.0:
          thermoplot(folder,"C"+"_"+str(i)+"_"+str(j),list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,ij]),yzero=0.0)
    thermoplot(folder,"B_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,43]),yzero=0.0)
    thermoplot(folder,"G_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,44]),yzero=0.0)
    thermoplot(folder,"E_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,45]),yzero=0.0)
    sys.exit()


  thermo = np.loadtxt(thermofile, comments="#", dtype=np.float)
  thermo[np.isnan(thermo)] = 0.0
  for i,cp in enumerate(thermo[:,6]):
    if cp > CpMax: break
    elif thermo[i,0] > tmelt: break

  thermo = thermo[0:i,:]
    
  Vstack=interpolate.splrep(thermo[:,0], thermo[:,1])
  V298 = float(interpolate.splev(T0, Vstack))
  Hstack=interpolate.splrep(thermo[:,0], thermo[:,4])
  H298 = float(interpolate.splev(T0, Hstack))
  threcord.update({"H298.15 (J/mol-atom)":round(H298,4)})
  Sstack=interpolate.splrep(thermo[:,0], thermo[:,3])
  S298 = float(interpolate.splev(T0, Sstack))
  threcord.update({"S298.15 (J/mol-atom/K)":round(S298,6)})

  zthermo.update({"temperature (K)":list(thermo[:,0])})
  zthermo.update({"atomic volume ($Angstrom^3$)":list(thermo[:,1])})
  thermoplot(folder,"atomic volume ($Angstrom^3$)",list(thermo[:,0]),list(thermo[:,1]))
  zthermo.update({"Gibbs energy (eV/atom)":list(thermo[:,2])})
  zthermo.update({"enthalpy (J/mol-atom)":list(thermo[:,4])})
  zthermo.update({"entropy (J/mol-atom/K)":list(thermo[:,3])})
  zthermo.update({"Cp (J/mol-atom/K)":list(thermo[:,6])})

  if fitCp:
    g,h,s,c,x=proStoichiometricCp()
  else:
    g,h,s,c,x=proStoichiometricG()

  threcord.update({"SGTE fitting":SGTErec})
  thermoplot(folder,"Gibbs energy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,2]*eVtoJ-H298),fitted=list(SGTE(x,g)), xT=list(x))
  thermoplot(folder,"enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), fitted=list(SGTEH(x,h)), xT=list(x))
  #thermoplot(folder,"enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), fitted=list(SGTE(x,g)+x*SGTES(x,s)), xT=list(x))
  thermoplot(folder,"entropy (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,3]),yzero=0.0, fitted=list(SGTES(x,s)), xT=list(x))

  zthermo.update({"LTC (1/K)":list(thermo[:,5])})
  thermoplot(folder,"LTC (1/K)",list(thermo[:,0]),list(thermo[:,5]),yzero=0.0)
  zthermo.update({"Cv (J/mol-atom/K)":list(thermo[:,14])})
  zthermo.update({"Cv,ion (J/mol-atom/K)":list(thermo[:,7])})
  Cele = [round(c,6) for c in thermo[:,14]-thermo[:,7]]
  zthermo.update({"Cele (J/mol-atom/K)":Cele})
  ncols = [6,14,7]
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]),fitted=list(SGTEC(x,c)), xT=list(x))
  ncols = [6,8]
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=300,expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=70,expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=100,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=1000,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=10000,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=300, expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=300, expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=70, expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=70, expt=expt, CoT=True)
  zthermo.update({"Debye temperature (K)":list(thermo[:,13])})
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0, xlim=70)
  zthermo.update({"bulk modulus (GPa)":list(thermo[:,15])})
  thermoplot(folder,"bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,15]),yzero=0.0)

  threcord.update({"zthermodynamic properies":zthermo})
  threcord.update({"Atomic volume at 298.15 K (Angstrom^3)":round(V298,6)})

  with open(vdos_e.replace('thermo/vdos_e','tplate/POSCAR'), 'r') as f:
    vvv = f.readlines()
  natom = sum([int(vv) for vv in vvv[6].split(' ') if vv!=""])
  structure.update({"number of atoms in POSCAR":natom})

  #natom = threcord.get("number of atoms in the primitive unit cell")

  with open(vdos_e.replace('thermo/vdos_e','thermo/data.in'), 'r') as f:
    vvv = f.readlines()
    Vfiles = []
    Pfiles = []
    volumes = []
    energies = []
    for vv in vvv[1:]:
      v = vv.split(' ')
      Pfiles.append("phonon/"+v[2].split('/')[-1].replace('\n', '').replace('"', ''))
      Vfiles.append(v[2].split('/')[-1].replace('\n', '').replace('"', ''))
      volumes.append(round(float(v[0])/natom,6))
      energies.append(round(float(v[1])/natom,6))
  structure.update({"Static vasp settings":Vfiles})
  structure.update({"phonon vasp settings and force constants":Pfiles})
  threcord.update({"volumes":volumes})
  threcord.update({"energies":energies})

  with open(dir0+'/E-V.dat','w') as f:
    for i,v in enumerate(volumes):
      f.write('{} {}\n'.format(v,energies[i]))
  #cmd = "YWfit -BMvol <"+dir0+'/E-V.dat | grep "f_expr(x) = "'
  ffun = "-Morse"
  cmd = "YWfit "+ffun+" <"+dir0+'/E-V.dat | grep "f_expr(x) = "'
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
  cwd = os.getcwd()
  os.chdir( dir0+"/figures")
  fitF = output.stdout
  with open('E-V.plt','w') as f:
    f.write('set terminal postscript landscape enhanced color "Times_Roman" 20\n')
    f.write('set encoding iso_8859_1\n')
    f.write('set pointsize 1.2\n')
    f.write('set size 0.95,0.95\n')
    f.write('set output "E-V.eps"\n')
    f.write('{}\n'.format(fitF))
    f.write('set key right bottom\n')
    f.write('set xlabel "atomic volume ($Angstrom^3$)\n')
    f.write('set ylabel "static energy (eV/atom)\n')
    f.write('plot "../E-V.dat" title "calculated" w p pt 7, \\\n')
    f.write('     f_expr(x) title "'+ffun+'" w l lt -1\n')
  #cmd = "gnuplot E-V.plt; convert -fuzz 100% -transparent white -rotate 90 -density 120x120 E-V.eps E-V.png"
  cmd = "gnuplot E-V.plt; convert -background white -alpha remove -rotate 90 -density 120x120 E-V.eps E-V.png"
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)

  figures.update({"static E-V curve": "figures/E-V.png"})
  threcord.update({"figures":figures})
  os.chdir(cwd)
  return Vfiles,Pfiles,g

def BMfitP(x,z):
  p = 0.0
  N = len(z)
  for n in range(N):
    p = p + float(N-n-1)*z[n]*x**(N-n+0.5)*(-2.0/3.0)
  return (-p)

def BMfit(v,p,g, T):
  global tPmax
  v = np.array(list(map(float,v)))
  p = np.array(list(map(float,p)))
  g = np.array(list(map(float,g)))
  f = g - p*v
  x = v**(-2.0/3.0)
  z = np.polyfit(x,f,BMvol)
  gf = np.poly1d(z)
  if (Debug==1):
    for i, vv in enumerate(v):
      print(vv, p[i], BMfitP(x[i],z), f[i], gf(x[i]))
   
  xx = x[0]
  xd = (x[len(x)-1] - x[0])*0.02
  pp = []
  gg = []
  vv = []
  for i in range(999):
    ppxx = BMfitP(xx,z)
    if (ppxx > tPmax*1.1):
      break
    pp.append(ppxx)
    vv.append(xx**(-1.5))
    gg.append(gf(xx)+ppxx*xx**(-1.5))
    xx = xx + xd
  try:
    s = interpolate.splrep(pp, gg)
    sv = interpolate.splrep(pp, vv)
  except ValueError:
    print("*******fetal ERROR: BMvol of order: ", BMvol, "  fetal fitting error at T= ", T)
    print(pp)
    print(gg)
    sys.exit()
  gx =interpolate.splev(txx, s)
  vx =interpolate.splev(txx, sv)
  if (Debug==1):
    for i, pp in enumerate(txx):
      print(pp*eVtoGPa, gx[i])
    for i, pp in enumerate(p):
      print(pp*eVtoGPa, g[i])
    sys.exit()
  class result:
    G = gx
    V = vx
  return(result)

def mkDict(line):
  rec = {}
  skiprec = False
  ss = str(line)[0:].replace("'","").split()
  ss = [k for k in ss if k != '']
  for nc,el in enumerate(ss):
    if isfloat(el):
      break
  if len(within)!=0:
    for el in ss[0:nc]:
      if el not in within:
        return True, None, None, None

  _sideal = 0
  _PN = ""
  i = nc*2+2
  while i < len(ss):
    #print ("ncx=", nc, ss[i])
    if ss[i] == "PQ":
      try:
        _PQ = float(ss[i+1])
        #threcord.update({"amount of imaginary phonon mode":float('{:.6f}'.format(_PQ))})
        Uncertainty.update({"amount of imaginary phonon mode":round(_PQ,6)})
        skiprec = _PQ >= PQ
        i += 1
        if skiprec:
          if not paper: print (ss[nc*2+1],"skipped, PQ=", ss[i])
          break
      except:
        skiprec = True
        print ("********Wrong record", ss)
        break
    elif ss[i] == "EQ":
      try:
        _EQ =  float(ss[i+1])
        #threcord.update({"0 K energy uncertainty (eV/atom)":float('{:.6f}'.format(_EQ))})
        Uncertainty.update({"0 K energy uncertainty (eV/atom)":round(_EQ,6)})
        skiprec = _EQ >= EQ
        i += 1
        if skiprec:
          print (ss[nc*2+1],"skipped, EQ=", ss[i])
          break
      except:
        skiprec = True
        print ("********Wrong record", ss)
        break
    elif ss[i] == "PN":
        _PN = ss[i+1].strip("/")
        threcord.update({"Phase name":_PN})
        i += 1
        mpid = ""
        try:
          mp = _PN.index("mp-")
          mpid = _PN[mp:].split('_')[0]
        except:
          pass
        structure.update({"mpid":mpid})
    elif ss[i] == "E0":
        #threcord.update({"static energy (eV/atom)":float('{:.6f}'.format(float(ss[i+1].strip("/"))))})
        threcord.update({"Static energy (eV/atom)":round(float(ss[i+1]),6)})
        i += 1
    elif ss[i] == "TT":
        Tup = float(ss[i+1])
        threcord.update({"Tmax":Tup})
        i += 1
        if Tup < Tupmax:
          print (ss[nc*2+1],"skipped, Tmax=", ss[i])
          skiprec = True
          break
    elif isfloat(ss[i]):
      _sideal = float(ss[i])
    i += 1

  if skiprec:
    return True, None, None, None

  if nc!=0:
    threcord.update({"Uncertainty":Uncertainty})
    space = ss[nc*2].strip("/").split("|")
    structure.update({"space group":int(space[0])})
    structure.update({"point group symmetry":space[1]})
    structure.update({"space group symmetry":space[2]})
    structure.update({"primitive unit cell formula":space[3]})
    elist, clist = formula2composition(space[3])
    pnatom = sum(clist)
    structure.update({"number of atoms in the primitive unit cell":int(pnatom)})
  
    tComponents = ss[0:nc]
    tnComponents = np.array(list(map(int,ss[nc:nc+nc])))
    natom = sum(tnComponents)
    tnComponents = tnComponents/natom
  
    Components = sorted(set(tComponents))
    nComponents = np.zeros(len(Components))
    for i0,el in enumerate(tComponents):
      ix = Components.index(el)
      nComponents[ix] = nComponents[ix] +  tnComponents[i0]
  
    compositions = []
    for i in range(len(Components)):
      compositions.append(int(0.1+natom*nComponents[i]))
    threcord.update({"Elements":Components})
    threcord.update({"Occupancies":list(compositions)})
    
  
    i = nc*2+2
    while i < len(ss):
      if ss[i] == "disordered":
        if i+1>=len(ss):
          _sideal = -sum(nComponents*np.log(nComponents))
        elif isfloat(ss[i+1]):
          i += 1
          if float(ss[i])<0.0:
            _sideal = -sum(nComponents*np.log(nComponents))
          else:
            _sideal = float(ss[i])
        else:
          _sideal = -sum(nComponents*np.log(nComponents))
      i += 1
    threcord.update({"Ideal mixing entropy (kB/atom)":_sideal})

  keys = threcord.keys()
  if nc==0: nc=-1
  vdos_e = str(ss[nc+nc+1]).replace('//','/')
  threcord.update({"Calculation date":str(datetime.datetime.fromtimestamp(os.path.getmtime(vdos_e)))})
  #threcord.update({"Calculation date":str(date.fromtimestamp(os.path.getatime(vdos_e)))})
  if _PN=="":
    try:
      _PN = [s for s in vdos_e.split('/') if s!=""][-3]
    except:
      _PN = "unknown"

  dir0 = _PN
  idx = 1
  while True:
    if not os.path.exists(dir0): break
    recordfile = dir0+"/record.json"
    newdir = False
    try:
      if os.path.exists(recordfile):
        with open(recordfile) as jsonfile:
          orec = json.load (jsonfile)
        okeys = orec.keys()
        for k in keys:
          v = threcord.get(k)
          for ok in okeys:
            if ok != k: continue
            newdir = orec.get(ok) != v
            if newdir: break 
          if k == "Static energy":
            if k in okeys:
              if abs(float(v-okeys.get(k))) < THR0: newdir = False
    except:
      pass

    if not newdir: break
    idx += 1
    dir0 = _PN+"#"+str(idx)

  oldPN = dir0
  newPN = PhaseName.get(_PN)
  if newPN != None: dir0 = newPN
  #print (_PN,dir0,PhaseName)

  threcord.update({"Phase name":dir0})

  """
  if paper:
    pname = threcord.get("primitive unit cell formula")
    n0 = 1
    for px in papers:
      if pname == px.split('#')[0]:
        n0 += 1
    if n0!=1:
      pname = pname +'#'+str(n0)
    papers.append(pname)
    dir0 = pname
  global start
  print ("thermo files extracting cost", time.time()-start)
  start = time.time()
  """
  return False, vdos_e, dir0, oldPN
  

def VASPResults(dir0,vdos_e,Vfiles, Pfiles, phdft="phonon"):
  hdir = vdos_e.replace('thermo/vdos_e','')
  
  natom = structure.get("number of atoms in POSCAR")
  pdir = vdos_e.replace('thermo/vdos_e',phdft)+'/'
  phdir = dir0+'/phonon'
  if not os.path.exists(phdir):
    os.mkdir(phdir)
  for ff in Vfiles:
    vdir = dir0+"/"+ff
    if not os.path.exists(vdir):
      os.mkdir(vdir)
    pvdir = phdir+"/"+ff
    if not os.path.exists(pvdir):
      os.mkdir(pvdir)
       
    vdos = hdir+"/"+ff+"/vdos.out"
    copyfile(vdos,vdir+'/vdos.out')
    print(vdos)

    poscar = hdir+"/"+ff+"/CONTCAR"
    if not os.path.exists(poscar):
      poscar = hdir+"/"+ff+"/Static.CON"
      if not os.path.exists(poscar):
        poscar = hdir+"/"+ff+"/POSCAR"
    copyfile(poscar,vdir+'/POSCAR')

    outcar = hdir+ff+"/OUTCAR"
    if os.path.exists(outcar):
      output = subprocess.run("grep POTCAR "+outcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)
    else:
      outcar = pdir+"/"+ff+"/OUTCAR.gz"
      if not os.path.exists(outcar):
        outcar = pdir+"/"+ff+"/Static.OUT.gz"
      output = subprocess.run("zgrep POTCAR "+outcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)

    #print (outcar)
    with open(vdir+'/POTCAR', "w") as text_file:
      text_file.write(POTCAR)

    doscar = hdir+"/"+ff+"/DOSCAR"
    if not os.path.exists(doscar):
      doscar = hdir+"/"+ff+"/Static.DOS.gz"
    ddoscar = doscar.split('/')[-1].replace("Static.DOS", "DOSCAR")
    copyfile(doscar, vdir+'/'+ddoscar)

    oszicar = hdir+"/"+ff+"/OSZICAR"
    if not os.path.exists(doscar):
      oszicar = hdir+"/"+ff+"/Static.OSZ"
    output = subprocess.run("grep E0= "+oszicar+" | tail -1 | awk '{print $5}'", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
    E0 = float(output.stdout)
    E0 /= natom
    with open(vdir+'/energy', "w") as text_file:
      text_file.write(str(E0))

    incars = fnmatch.filter(os.listdir(hdir+"/"+ff), 'INCAR_*')
    if len(incars)==0:
      INCAR = hdir+"/tplate/INCAR.Static"
    else:
      INCAR = hdir+"/"+ff+'/'+incars[0]
      d0 = os.path.getmtime(INCAR)
      for ii in range(1,len(incars)):
        d1 = os.path.getmtime(hdir+"/"+ff+'/'+incars[ii])
        if d1 > d0:
          INCAR = hdir+"/"+ff+'/'+incars[ii]
          d1 = d0
    copyfile(INCAR, vdir+'/INCAR')
    copyfile(hdir+"/tplate/KPOINTS", vdir+'/KPOINTS')

    sposcar = pdir+ff+"/POSCAR"
    copyfile(sposcar, pvdir+'/POSCAR')
    soutcar = pdir+ff+"/OUTCAR"
    if os.path.exists(soutcar):
      output = subprocess.run("grep POTCAR "+soutcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)
    else:
      soutcar = pdir+"/"+ff+"/OUTCAR.gz"
      output = subprocess.run("zgrep POTCAR "+soutcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)

    with open(pvdir+'/POTCAR', "w") as text_file:
      text_file.write(POTCAR)

    incars = fnmatch.filter(os.listdir(pdir+"/"+ff), 'INCAR_*')
    if len(incars)==0:
      INCAR = pdir+"/tplate/INCAR"
    else:
      INCAR = pdir+"/"+ff+'/'+incars[0]
      d0 = os.path.getmtime(INCAR)
      for ii in range(1,len(incars)):
        d1 = os.path.getmtime(pdir+"/"+ff+'/'+incars[ii])
        if d1 > d0:
          INCAR = pdir+"/"+ff+'/'+incars[ii]
          d1 = d0
    copyfile(INCAR, pvdir+'/INCAR')
    copyfile(pdir+"/tplate/KPOINTS", pvdir+'/KPOINTS')

    sposcar = pdir+ff+"/POSCAR"
    sxml = pdir+ff+"/vasprun.xml"
    if not os.path.exists(sxml):
      sxml = pdir+ff+"/vasprun.xml.gz"

    cwd = os.getcwd()
    os.chdir( pvdir )
    cmd = 'vasp_fij -outc '+soutcar+" -xml "+sxml+" -conc "+sposcar + " >& /dev/null"
    os.system(cmd)
    os.chdir( cwd )

  global start
  print ( round(time.time()-start,3), "Secs. costed in VASP files extracting")
  start = time.time()

def extractGph():
  phononmode = {}
  with open("symmetry.out", "r") as f:
    lines = f.readlines()
  i = 0
  while i < len(lines):
    ss = [s for s in lines[i].strip().split(' ') if s!='']
    if len(ss) >= 5:
      if ss[2] == "Modes" and ss[4] in ["silent_mode", "raman_active", "ir_active"]:
        mode = []
        for ii in range(int(ss[0])):
          i += 1
          mm = [s for s in lines[i].strip().replace('(',' ').replace(')',' ').split(' ') if s!='']
          mode.append(float(mm[3]))
        phononmode.update({ss[1]+" ( "+ss[4]+" )": sorted(mode)})
    i += 1
  threcord.update({"gamma point phonons (cm-1) ":phononmode})
    
def Phonon298(dir0, pvdos=False):
  V298 = threcord.get("Atomic volume at 298.15 K (Angstrom^3)")
  phdir298 = dir0 + '/phonon298.15K'
  if not os.path.exists(phdir298):
    os.mkdir(phdir298)
  volumes = threcord.get("volumes")
  i1 = 0
  for ii,vv in enumerate(volumes):
    if float(vv) < V298:
      i1 += 1
  i1 -= 1
  i1 = max(i1, 0)
  i1 = min(i1, len(volumes)-2)
  dV = float(volumes[i1+1]) - float(volumes[i1])
  ff1 = (float(volumes[i1+1]) - V298)/dV
  cmd = "Ymix -f "+str(ff1)+ " " + dir0+'/'+Pfiles[i1]+"/superfij.out " + " " + dir0+'/'+Pfiles[i1+1]+"/superfij.out >"+phdir298+"/superfij.out"
  print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)

  cwd = os.getcwd()
  os.chdir( phdir298 )

  cmd = "Yphon -tranI 2 -eps -nqwave "+ str(nqwave)+ " <superfij.out"
  print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  cmd = "gnuplot vdos.plt; convert background white -alpha remove -rotate 90 -density 120x120 vdos.eps vdos.png"
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  figures = threcord.get("figures")
  figures.update({"phonon DOS at 298.15 K": "phonon298.15K/vdos.png"})

  if pvdos:
    cmd = "Yphon -tranI 2 -eps -pvdos -nqwave "+ str(nqwave/4)+ " <superfij.out"
    print(cmd)
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    cmd = "gnuplot pvdos.plt; convert background white -alpha remove -rotate 90 -density 120x120 pvdos.eps pvdos.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    figures = threcord.get("figures")
    figures.update({"generalized phonon DOS at 298.15 K": "phonon298.15K/pvdos.png"})
  os.chdir( cwd )

  #ngroup = threcord.get("structure").get("space group")
  ngroup = structure.get("space group")
  dfile = ""
  if ngroup>=1 and ngroup<=2:
    dfile = home+"/bin/hbin/pycode/data/dfile.tri"
    structure.update({"crystal system": "Triclinic"})
  elif ngroup>=3 and ngroup<=15:
    dfile = home+"/bin/hbin/pycode/data/dfile.mon"
    structure.update({"crystal system": "Monoclinic"})
  elif ngroup>=16 and ngroup<=74:
    dfile = home+"/bin/hbin/pycode/data/dfile.oth"
    structure.update({"crystal system": "Orthorhombic"})
  elif ngroup>=75 and ngroup<=142:
    dfile = home+"/bin/hbin/pycode/data/dfile.tet"
    structure.update({"crystal system": "Tetragonal"})
  elif ngroup>=143 and ngroup<=167:
    dfile = home+"/bin/hbin/pycode/data/dfile.rho"
    structure.update({"crystal system": "Trigonal"})
  elif ngroup>=168 and ngroup<=194:
    dfile = home+"/bin/hbin/pycode/data/dfile.hcp"
    structure.update({"crystal system": "Hexagonal"})
  elif ngroup>=195 and ngroup<=220:
    dfile = home+"/bin/hbin/pycode/data/dfile.scc"
    structure.update({"crystal system": "Cubic"})
  elif ngroup>=221 and ngroup<=224:
    dfile = home+"/bin/hbin/pycode/data/dfile.bcc"
    structure.update({"crystal system": "Cubic({bcc})"})
  elif ngroup>=225 and ngroup<=230:
    dfile = home+"/bin/hbin/pycode/data/dfile.fcc"
    structure.update({"crystal system": "Cubic({fcc})"})

  if dfile != "":
    dfile0 = dfile.split('/')[-1]
    copyfile(dfile,phdir298+'/'+dfile0)
    cwd = os.getcwd()
    os.chdir( phdir298 )
    cmd = 'timeout 6 pos2s Symmetry.pos -THR 1.e-4 >&symmetry.out'
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)

    Gph = os.path.exists("symmetry.mode")
    if Gph:
      cmd = "Yphon -Gfile symmetry.mode -tranI 2 -eps -pdis "+dfile0+ " <superfij.out >symmetry.out"
    else:
      cmd = "Yphon -tranI 2 -eps -pdis "+dfile0+ " <superfij.out >symmetry.out"

    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    if Gph:
      extractGph()

    cmd = "gnuplot vdis.plt; convert -background white -alpha remove -rotate 90 -density 120x120 vdis.eps vdis.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    figures.update({"phonon dispersion at 298.15 K": "phonon298.15K/vdis.png"})
    os.chdir( cwd )

  threcord.update({"figures":figures})
  
  cwd = os.getcwd()
  os.chdir( phdir298 )
  cfile = ["findsym.log","run.log","exactQ.out","run","symmetry.out"]
  for f in cfile:
    if os.path.exists(f):
      os.remove(f)
  os.chdir( cwd )
    
  global start
  print (round(time.time()-start,3), "Secs. costed in calculations of phonon properties at 298.15K")
  start = time.time()

#outs = ["space group", "point group symmetry", "space group symmetry"]
outs = ["space group","space group symmetry"]
def addpapers(g,formula,pname):
  g[0] += float(threcord.get("H298.15 (J/mol-atom)"))
  if threcord.get("Ideal mixing entropy (kB/atom)")!=None:
    g[1] -= R*float(threcord.get("Ideal mixing entropy (kB/atom)"))

  sys.stdout.write("{},{}".format(pname,formula))
  for ss in outs:
    sys.stdout.write(",{}".format(structure.get(ss)))
  sys.stdout.write(",{:.6g}".format(g[0]))
  sys.stdout.write(",{:.6g}".format(g[1]))
  sys.stdout.write(",{:.5g}".format(g[2]))
  sys.stdout.write(",{:.3e}".format(g[3]))
  sys.stdout.write(",{:.3e}".format(g[4]))
  sys.stdout.write(",{:.3e}".format(g[5]))

  try:
      pq = 100.*threcord.get("Uncertainty").get("amount of imaginary phonon mode")
      eq = eVtoJ*threcord.get("Uncertainty").get("0 K energy uncertainty (eV/atom)")
  except:
      pq = 0.0
      eq = 0.0
  warning = SGTErec.get("G-H298.15 (J/mol-atom)")[1].get("error")
  #print(pq,eq,warning)
  #sys.stdout.write(",{:.1f},{:.0f},{:.0f},".format(pq,eq,warning))
  sys.stdout.write(",{:.1f},".format(pq))

  if warning > 99:
    sys.stdout.write("   **********WARNING, fittit error is too large! {}".format(warning)) 
  if abs(g[2]) > 32:
    sys.stdout.write("   **********WARNING, Cp at room condition is abnormal! {}".format(abs(g[2]))) 
  sys.stdout.write("\n")


markers=['o', 'v', 'd', '^', '<', '>', 's', '*', 'x', '+', '1', '2']

home = "/global/u2/y/yiwang62/"
k_B = 8.6173303e-5
R = 8.3144598

eVtoGPa = 160.21766208 
eVtoJ = 96486.9
THRE0 = 1.e-5
nqwave = 2.e6
nqwave = 1.e6

T0 = 298.15
update = True

input_within = False #""" key to cotrol the within input"""
formula_within = "" #"""chemical formula"""
within = []

PQ = 0.075
EQ = 0.015
PQ = 0.01
EQ = 0.01
CpMax = 50.
Tupmax = 2000.0
start = time.time()
threcord = {}
figures = {}
zthermo = {}
SGTErec = {}
structure = {}
expt = None
Uncertainty = {}
debug = False
fitCp = True
paper = True
phdft = "phonon"
justplot = None
pvdos = False

phases = []
papers = []
PhaseName = {}

def plotAPI(thermofile, volumes, energies, expt=None, xlim=None):
  phasedir = [substr for substr in thermofile.split('/') if substr!=""]
  phasedir = ('/').join(phasedir[0:-1])
  if phasedir=="": phasedir="."
  folder = phasedir+"/figures/"
  print("All figures have been outputed into: ", folder, "  with T uplimt:", xlim, "\n\nEnjoy!\n")
  if not os.path.exists(folder):
    os.mkdir(folder)
  thermoplot(folder,"0 K total energies (eV/atom)",volumes, energies)

  thermo = np.loadtxt(thermofile, comments="#", dtype=np.float)
  thermo[np.isnan(thermo)] = 0.0
  for i,cp in enumerate(thermo[:,6]):
    if cp > CpMax: 
      thermo = thermo[0:i,:]
      break
    
  """
  f2=interpolate.splrep(thermo[:,0], thermo[:,1])
  V298 = float(interpolate.splev(T0, Vstack))
  Hstack=interpolate.splrep(thermo[:,0], thermo[:,4])
  H298 = float(interpolate.splev(T0, Hstack))
  Sstack=interpolate.splrep(thermo[:,0], thermo[:,3])
  S298 = float(interpolate.splev(T0, Sstack))
  """
  f2=interp1d(thermo[:,0], thermo[:,1])
  V298 = f2(T0)
  f2=interp1d(thermo[:,0], thermo[:,4])
  H298 = f2(T0)
  f2=interp1d(thermo[:,0], thermo[:,3])
  S298 = f2(T0)

  threcord.update({"H298.15 (J/mol-atom)":H298})
  threcord.update({"S298.15 (J/mol-atom/K)":S298})

  zthermo.update({"temperature (K)":list(thermo[:,0])})
  zthermo.update({"atomic volume ($Angstrom^3$)":list(thermo[:,1])})
  zthermo.update({"Gibbs energy (eV/atom)":list(thermo[:,2])})
  zthermo.update({"enthalpy (J/mol-atom)":list(thermo[:,4])})
  zthermo.update({"entropy (J/mol-atom/K)":list(thermo[:,3])})
  zthermo.update({"Cp (J/mol-atom/K)":list(thermo[:,6])})
  proStoichiometricCp()
  with open(folder + '/../record.json', 'w') as fp:
    myjsonout(SGTErec, fp, indent="", comma="")
  myjsonout(SGTErec, sys.stdout, indent="", comma="")

  thermoplot(folder,"atomic volume ($Angstrom^3$)",list(thermo[:,0]),list(thermo[:,1]), xlim=xlim)
  thermoplot(folder,"Gibbs energy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,2]*eVtoJ-H298), xlim=xlim)
  thermoplot(folder,"enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), xlim=xlim)
  thermoplot(folder,"entropy (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,3]),yzero=0.0, xlim=xlim)

  thermoplot(folder,"LTC (1/K)",list(thermo[:,0]),list(thermo[:,5]),yzero=0.0, xlim=xlim)
  ncols = [6,8]
  thermoplot(folder,"heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), expt=expt, xlim=xlim)
  thermoplot(folder,"heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=300,expt=expt)
  thermoplot(folder,"heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=100,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), elonly=300, expt=expt, CoT=True)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,10]),yzero=0.0, xlim=xlim)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,10]),yzero=0.0, xlim=70)
  thermoplot(folder,"bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,9]),yzero=0.0,xlim=xlim)
  thermoplot(folder,"Seebeck coefficients (μV/K)",list(thermo[:,0]),list(thermo[:,16]),xlim=xlim)
  thermoplot(folder,"Lorenz number ($WΩK^{−2}$)",list(thermo[:,0]),list(thermo[:,17]),xlim=xlim)


if __name__ == '__main__':
    count = 1
    while (count < len(sys.argv)):
      if (input_within):
        if sys.argv[count].startswith('-'):
          input_within = False
        else:
          formula_within = formula_within+sys.argv[count]
          count = count + 1
          if (count > len(sys.argv)):
            break
          continue
      if (sys.argv[count] == "-pvdos"):
        pvdos = True
      if (sys.argv[count] == "-within"):
        input_within = True
      elif (sys.argv[count] == "-expt"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        with open(sys.argv[count], encoding='utf-8') as element_json_file:
          expt = json.load(element_json_file)
      elif (sys.argv[count] == "-phdft"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        phdft = sys.argv[count]
      elif (sys.argv[count] == "-T0"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        T0 = float(sys.argv[count])
      elif (sys.argv[count] == "-phasename"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        phasename = str(sys.argv[count])
      elif (sys.argv[count] == "-cpmax"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        CpMax = float(sys.argv[count])
      elif (sys.argv[count] == "-THRE0"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        THRE0 = float(sys.argv[count])
      elif (sys.argv[count] == "-PQ"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        PQ = float(sys.argv[count])
      elif (sys.argv[count] == "-EQ"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        EQ = float(sys.argv[count])
      elif (sys.argv[count] == "-Tupmax"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        Tupmax = float(sys.argv[count])
      elif (sys.argv[count] == "-nqwave"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        nqwave = float(sys.argv[count])
      elif (sys.argv[count] == "-debug"):
        debug = True
      elif (os.path.exists(sys.argv[count])):
        justplot=sys.argv[count]
      else:
        print ("*******Unknown option", sys.argv[count])
      count = count + 1
    
    if formula_within!="":
      within = formula2elist(formula_within)
      print ("data to be extracted within ",within)
    
    #print (phasename)
    #if True:
    try:
      with open (phasename,'r') as f:
        lines = f.readlines()
      for ll in lines:
        line = ll.strip('\n').replace(',', ' ').replace(':', ' ')
        ss = [s.strip() for s in line.split(' ') if s.strip()!='']
        if len(ss)>1: PhaseName.update({ss[0]:ss[1]})
    except:
      pass
    """
    """
    
    
    if justplot==None: lines = sys.stdin.readlines()
    else: lines = [justplot]
    
    sys.stdout.write("G(T)=a+b*T+c*T*Ln(T)+d*T*T+e*T*T*T+f/T (J/mol-atom)\n")
    sys.stdout.write("Phase,comp")
    for ss in outs:
      sys.stdout.write(",{}".format(ss))
    #sys.stdout.write(",a,b,c,d,e,f,PQ(%),EQ(J),GQ(J),\n".format(ss))
    sys.stdout.write(",a,b,c,d,e,f,PQ(%),\n".format(ss))
    
    for line in lines:
      if line.strip()=="": continue
      skip,vdos_e,dir0,oldPN = mkDict(line)
      if skip: continue
      #phases.append([threcord.get("phase name"),threcord.get("mpid")])
      phases.append(threcord.get("Phase name"))
      nphases += 1
      #print (line.strip())
    
      if not os.path.exists(dir0):
        os.mkdir(dir0)
      else:
        if not update: continue
    
      Vfiles,Pfiles,g = Genergy(vdos_e,dir0)
      addpapers(g,dir0,oldPN)
    
      if not debug:
        VASPResults(dir0,vdos_e,Vfiles, Pfiles, phdft=phdft)
        try:
          Phonon298(dir0, pvdos=pvdos)
        except:
          pass
      threcord.update({"structure":structure})
      #threcord.delete({"number of atoms in the primitive unit cell")
    
      with open(dir0 + '/record.json', 'w') as fp:
        myjsonout(threcord, fp, indent="", comma="")
    
    print ("\n", phases, "\n")
    print ("\n", nphases, "phases extracted\n")
