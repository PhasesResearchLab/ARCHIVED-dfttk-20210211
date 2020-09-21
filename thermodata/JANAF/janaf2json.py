import sys
from dfttk.analysis.ywutils import reduced_formula, formula2composition

def endoffile(lines):
    form = None
    for l in range(len(lines)):
        ss = [f for f in lines[l].strip().split('\t') if f!='']
        if len(ss) < 1: continue
        formula = ss[-1].split('(')
        if len(formula) < 2: continue
        try:
            els, natom = formula2composition(formula[0], False)
            form = reduced_formula(formula[0])
            lines = lines[l+1:]
            break
        except:
            continue
    if form==None: return True
    else: return False


def endofrec(line0,line1):
    ss = [f for f in line0.strip().replace('\t', ' ').split(' ') if f!='']
    try:
        float(ss[0])
    except:
        sys.stdout.write('\n        ]\n    }')
        return True
    #print(ss)
    if 'CRYSTAL' in ss:
        sys.stdout.write('\n        ]\n    }')
        return True
    elif 'LIQUID' in ss:
        sys.stdout.write('\n        ]\n    }')
        return True
    try:
        if len(line1.strip())==0:
            sys.stdout.write('\n        ]\n    }')
            return True
    except:
        sys.stdout.write('\n        ]\n    }')
        return True

    ss = [f for f in line1.strip().split('\t') if f!='']
    if len(ss)==0:
        sys.stdout.write('\n        ]\n    }')
        return True

    try:
        float(ss[0])
    except:
        sys.stdout.write('\n        ]\n    }')
        return True

def gatdata(lines):
    for k in range(1, len(lines)):
        ss = [f for f in lines[k].strip().split('\t') if f!='']
        if len(ss)==0: continue
        sys.stdout.write('{}, {}'.format(float(ss[0]), float(ss[1])))
        if endofrec(lines[k], lines[k+1]):
            """
            try:
                if lines[k+1].strip()=='': sys.stdout.write(',\n')
            except:
                sys.stdout.write('\n]')
            """
            kk = k
            break
        else:
            sys.stdout.write(',\n')
    for k in range(kk+1, len(lines)):
        try:
            if lines[k].strip()!='': 
                if not endoffile(lines[k:]):
                    sys.stdout.write(',\n')
                return k
        except:
            sys.stdout.write('\n]')
            return k
    return kk+1
        

lines = sys.stdin.readlines()

#l = 0
#while l < len(lines):
sys.stdout.write('[')
while True:
    form = None
    for l in range(len(lines)):
        ss = [f for f in lines[l].strip().split('\t') if f!='']
        if len(ss) < 1: continue
        formula = ss[-1].split('(')
        if len(formula) < 2: continue
        try:
            els, natom = formula2composition(formula[0], False)
            form = reduced_formula(formula[0])
            lines = lines[l+1:]
            break
        except:
            continue
    if form==None: break
    #print("eeeeeeeee", formula)

    sys.stdout.write('\n    {\n        "Author":"Chase(JANAF), Malcolm W. Chase, Jr. NIST-JANAF Thermochemical Tables. Washington, DC : New York :American Chemical Society ; American Institute of Physics for the National Institute of Standards and Technology, 1998.",\n')
    sys.stdout.write('        "Compound": "{}",\n'.format(form))
    sys.stdout.write('        "Unit": "J/K",\n')
    sys.stdout.write('        "natom": {},\n'.format(int(sum(natom))))
    sys.stdout.write('        "property": "heat capacity",\n')
    sys.stdout.write('        "data": [\n')

    try:
        kk = gatdata(lines)
        lines = lines[kk:]
        done = True
        for line in lines:
            if line.strip()!='': 
                done = False
                break
        if done: break
    except:
        break
sys.stdout.write('\n]\n')
