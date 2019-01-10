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

pdb_inp = iotbx.pdb.input(file_name="acid060.pdb")
hierarchy = pdb_inp.construct_hierarchy()

def cutf(R):
    if angular_cutoff >= R:
       Fc = 0.5*math.cos(math.pi*R/radial_cutoff)+0.5
    else:
       Fc = 0
    return Fc

def Rpart(Rj, Rs):
    GmRs = []
    n = 4.0
    for b in Rs:
        GmR = 0
        for a in Rj:
            f = cutf(a)
            if a >= 0:
                GmR = math.exp(- n * ((a - b) ** 2)) * f
        GmRs.append(GmR)
    return GmRs

def Apart(zeta, Rj, Rk):
    l = 8.00
    zetas = ts_values
    n = 4.0
    Rs = angular_rs_values
    GmAs = []
    for a in Rs:
        for b in zetas:
            GmA = 0
            for c in Rj:
                for d in Rk:
                    if c > 0:
                        if d > 0:
                            fj = cutf(c)
                            fk = cutf(d)
                            for e in zeta:
                                if e == 0:
                                    GmA = 0
                                else:
                                    GmA = GmA + ((( 1+math.cos(e - b)) )** l) * \
                                          math.exp(- n *( (((d + c) / 2) - a) ** 2)) * fj * fk
            GmA = GmA * (2 ** (1-l))
            GmAs.append(GmA)
    return GmAs
def Read(atom_types):
    GmRs = []
    target = open('text4.txt', 'w')
    for a1 in hierarchy.atoms():
        e1 = a1.element.upper()
        if e1 in atom_types:
            X = e1
            GmRs.append(X)
            target.write(str(X))
            target.write('\n')
            for J in atom_types:
                Rj = []
                for a2 in hierarchy.atoms():
                    e2 = a2.element.upper()
                    if e2 == J:
                        R = a1.distance(a2)
                        Rj.append(R)
                        for K in atom_types:
                            Rk = []
                            zeta = []
                            for a3 in hierarchy.atoms():
                                e3 = a3.element.upper()
                                if e3 == K:
                                    R = a1.distance(a2)
                                    Rk.append(R)
                                    ts = a1.angle(a2, a3, deg = True)
                                    if ts == None:
                                        ts = 0.0
                                    zeta.append(ts)
                            GmR = Apart(zeta, Rj, Rk)
                            target.write(str(J))
                            target.write(str(K))
                            target.write(str(GmR))
                            target.write('\n')
                            GmRs.append(J)
                            GmRs.append(K)
                            GmRs.append(GmR)
    target.close()
    print("Bingo", "*" * 50 )
    return GmRs

Read(atom_type)

