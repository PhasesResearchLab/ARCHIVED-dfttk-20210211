import h5py
import numpy as np
import sys
#filename = "/gpfs/scratch/yyx5048/CRO_data/phononpy/force_constants.hdf5"
filename = "force_constants.hdf5"

def get_alist():
  with h5py.File(filename, "r") as f:
    a_group_key = list(f.keys())[1]
    data = list(f[a_group_key])
  return data


def second_round(alist, pcell, scell):
  with h5py.File(filename, "r") as f:
    a_group_key = list(f.keys())[0]
    data = f[a_group_key]
    natom = data.shape[0]
    natoms = data.shape[1]
    #print(natom,natoms)
    hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)
    #print(data.shape, scell.shape)
    for cc in range(scell.shape[1]):
      for ix in range(len(pcell)):
        ii = scell[ix,cc]
        #print (ii,ix,cc)
        for jj in range(natoms):
            jx = None
            for bb in range(scell.shape[1]):
                for z in range(len(pcell)):
                    if jj==scell[z,bb]:
                        bx = bb-cc
                        if bx <0: bx += scell.shape[1]
                        jx = scell[z,bx]
                        break
                    if jx is not None: break
                if jx is not None: break
            for x in range(3):
                for y in range(3):
                    try:
                        hessian_matrix[ii*3+x, jj*3+y] = -data[ix,jx,x,y]
                    except:
                        print ("OUt of 2222",ii, jj, ix, jx)
                        sys.exit()
    #sys.exit()
    for xx in range(natoms*3):
        for yy in range(natoms*3-1):
            sys.stdout.write('{} '.format(hessian_matrix[xx,yy]))
        sys.stdout.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))


def first_round():
  with h5py.File(filename, "r") as f:
    a_group_key = list(f.keys())[0]
    data = f[a_group_key]
    natom = data.shape[0]
    natoms = data.shape[1]
    #print(natom,natoms)
    hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)
    for ii in range(natoms):
        for jj in range(natoms):
            for x in range(3):
                for y in range(3):
                    if ii < natom:
                        hessian_matrix[ii*3+x, jj*3+y] = -data[ii,jj,x,y]
                    else:
                        idx = ii%natom
                        jshift = ii//natom
                        jdx = (jj + jshift)%natoms
                        #print (idx, jdx)
                        hessian_matrix[ii*3+x, jj*3+y] = -data[idx,jdx,x,y]
    for xx in range(natoms*3):
        for yy in range(natoms*3-1):
            sys.stdout.write('{} '.format(hessian_matrix[xx,yy]))
        sys.stdout.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))


def get_natoms():
  with h5py.File(filename, "r") as f:
    a_group_key = list(f.keys())[0]
    data = f[a_group_key]
    natom = data.shape[0]
    natoms = data.shape[1]
  print (natom, natoms)


count = 1
amap = ""
while (count < len(sys.argv)):
  if (sys.argv[count] == "-map"):
    count = count + 1
    if (count > len(sys.argv)):
      break
    amap = sys.argv[count]
  elif (sys.argv[count] == "-getnatoms"):
    get_natoms()
    sys.exit()
  else:
    filename = sys.argv[count]
  count = count + 1

if amap != "":
  with open(amap, "r") as f:
    lines = f.readlines()
  els = []
  line = lines[0]
  field = [f for f in line.strip().split(' ') if f!=""]
  pcell = []
  scell = np.zeros((len(lines), len(field)-3), dtype=int)
  for i,line in enumerate(lines):
    field = [f for f in line.strip().split(' ') if f!=""]
    pcell.append(field[1])
    els.append(field[2])
    for j in range(0,len(field)-3):
      scell[i,j] = field[j+3]
  alist = get_alist()
  if len(alist)!=len(pcell):
    print ("ERROR: incorrect primitive unit cell")
    sys.exit()
  for s in alist:
    if s not in scell[:,0]:
      print ("ERROR: alist not covered by scell[:,0]")
      sys.exit()
  for s in scell[:,0]:
    if s not in alist:
      print ("ERROR: scell[:,0] not covered by alist")
      sys.exit()

  second_round(alist, pcell, scell)
  sys.exit()

first_round()
