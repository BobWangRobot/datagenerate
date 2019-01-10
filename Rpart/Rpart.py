import iotbx.pdb
import numpy
import math

rs_values = [ 0.900000,1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
angular_rs_values = [0.900000, 2.200000]
atom_type = [' C', ' H',' O']
radial_cutoff = 5.2
angular_cutoff = 3.5
radial_nu = 8
radial_angular_nu = 8
angular_zeta = 8

pdb_inp = iotbx.pdb.input(file_name="acid300.pdb")
hierarchy = pdb_inp.construct_hierarchy()

def cutf(R):
    if radial_cutoff >= R:
       Fc = 0.5*math.cos(math.pi*R/radial_cutoff)+0.5
    else:
       Fc = 0
    return Fc

def Rpart(Rj):
    GmRs = []
    n = 4.0
    Rs = rs_values
    for b in Rs:
        GmR = 0
        for a in Rj:
            f = cutf(a)
            if a > 0:
                GmR = math.exp(- n * ((a - b) ** 2)) * f
        GmRs.append(GmR)
    return GmRs

def Read(atom_types):
    GmRs = []
    target = open('text1.txt', 'w')
    for a1 in hierarchy.atoms():
        e1 = a1.element.upper()
        if e1 in atom_types:
            X = e1
            GmRs.append(X)
            target.write(str(X))
            target.write('\n')
            for F in atom_types:
                Rj = []
                for a2 in hierarchy.atoms():
                    e2 = a2.element.upper()
                    if e2 == F:
                        R = a1.distance(a2)
                        Rj.append(R)
                GmR = Rpart(Rj)
                target.write(str(F))
                target.write(str(GmR))
                target.write('\n')
                GmRs.append(F)
                GmRs.append(GmR)
    target.close()
    print("Bingo", "*" * 50 )
    return GmRs

Read(atom_type)

