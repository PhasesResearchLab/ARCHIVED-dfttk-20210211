#python thermodata.py */*.json
import json
import math
import sys
import numpy as np
from scipy.interpolate import interp1d
from dfttk.analysis.ywplot import Myjsonout


def angle(i,T,C,fac=1.0):
    vector_1 = [T[i] - T[i-1], fac*(C[i] - C[i-1])]
    vector_2 = [T[i+1] - T[i], fac*(C[i+1] - C[i])]
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = min(np.dot(unit_vector_1, unit_vector_2),1.0)
    angle = math.degrees(np.arccos(dot_product))
    #print (vector_1,vector_2,dot_product,angle)
    return angle


def average_angle(T,C):
    fac = max(T)/max(C)
    tot = 0.
    n = 0
    for i in range(1,len(T)-2):
        tot += angle(i,T,C,fac)
        n += 1
    return tot/n, fac


def newset(t,c, Tmax):
    dT = Tmax/30
    dT = min(dT,90)
    nT = int((t[-1]-t[0])/dT+.5)
    if nT <3: nT = 3
    tn = np.linspace(t[0], t[-1], nT)
    f2 = interp1d(t,c)
    cn = f2(tn)
    #print(t[0],t[-1])
    #print(tn[0],tn[-1])
    return tn,cn
    

def unique(T,C):
    _T = []
    _C = []
    for i, x in enumerate(T):
        if x in _T: continue
        _T.append(T[i])
        _C.append(C[i])
    return _T,_C


def reform(rec):
    _T,_C = unique(np.array(rec[0::2]), np.array(rec[1::2]))
    thr,fac = average_angle(_T,_C)
    #print ("thr=",thr,fac)
    Tmax = max(_T)
    T = []
    C = []
    t=[_T[0]]
    c=[_C[0]]
    while True:
        Done = True
        for i in range(1,len(_T)-2):
            t.append(_T[i])
            c.append(_C[i])
            if (_C[i]-_C[i-1])*(_C[i+1]-_C[i])<0 or \
                angle(i,_T,_C,fac)>30*thr:
                #print (_T[i-1], _T[i], _T[i+1])
                #print (_C[i-1], _C[i], _C[i+1])
                newT, newC = newset(t,c,Tmax)
                T.extend(newT)
                C.extend(newC)
                _T = _T[i+1:]
                _C = _C[i+1:]
                t=[_T[0]]
                c=[_C[0]]
                Done = False
                break
        if Done:
            t.append(_T[-1])
            c.append(_C[-1])
            newT, newC = newset(t,c, Tmax)
            T.extend(newT)
            C.extend(newC)
            break

    data = np.zeros((len(T)*2), dtype=float)
    data[0::2] = T
    data[1::2] = C
    return list(data)


record = []
for jj in range (1,len(sys.argv)):
    print (sys.argv[jj])
    with open (sys.argv[jj], "r") as fp:
        orec = json.load (fp)
        record.extend(orec)

for i,r in enumerate(record):
    if r['Author'].startswith("Andersson(CALPHAD)"):
        record[i]['data'] = reform(record[i]['data'])
        #for j in range(0, len(record[i]['data']),2):
        #    print(record[i]['data'][j], record[i]['data'][j+1])

    
with open ("ExptData.json", "w") as out:
    Myjsonout(record,out)

print ("\n", len(record), "records handled.\n")

with open ("ExptData.json", "r") as fp:
    orec = json.load (fp)

