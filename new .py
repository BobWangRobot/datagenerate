import iotbx.pdb
import numpy
import math

rs_values = [ 0.900000,1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
angular_rs_values = [0.900000, 2.200000]
radial_cutoff = 5.2
angular_cutoff = 3.5
radial_nu = 8
radial_angular_nu = 8
angular_zeta = 8



def Atome_classify(hierarchy):
    atom_elements = {}
    for a in hierarchy.atoms():
        e = a.element.upper().strip()
        atom_elements.setdefault(e, [])
        atom_elements[e].append(a)
    return atom_elements

def cutf(R):
    if radial_cutoff >= R:
       Fc = 0.5*math.cos(math.pi*R/radial_cutoff)+0.5
    else:
       Fc = 0
    return Fc

def Req(Rj):
    GmRs = []
    n = 4.0
    Rs = rs_values
    for b in Rs:
        GmR = 0
        for a in Rj:
            f = cutf(a)
            GmR += math.exp(- n * ((a - b) ** 2)) * f
        GmRs.append(GmR)
    return GmRs

def Rpart(hierarchy, file_name):
    print len(hierarchy.atoms())
    assert 0
    atom_types = get_atom_elements(hierarchy)
    print atom_types
    Fake = Atome_classify(hierarchy)
    target = open(file_name, 'w')
    for a1 in hierarchy.atoms():
        GmRs = []
        e1 = a1.element.upper().strip()
        assert e1 in atom_types, '%s not in %s' % (e1, atom_types)
        X = e1
        target.write(str(X))
        target.write('\n')
        for e2, atoms2 in Fake.items():
            Rj = []
            print e2
            print atoms2
            for a2 in atoms2:
                R = a1.distance(a2)
                Rj.append(R)
            GmR = Req(Rj)
            target.write('H')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
        assert 0
        Rj = []
        for a3 in C:
            R = a1.distance(a3)
            Rj.append(R)
        GmR = Req(Rj)
        target.write('C')
        target.write(str(GmR))
        target.write('\n')
        GmRs.append(GmR)
        Rj = []
        for a4 in O:
            R = a1.distance(a4)
            Rj.append(R)
        GmR = Req(Rj)
        target.write('O')
        target.write(str(GmR))
        target.write('\n')
        GmRs.append(GmR)
        target.write(str(GmRs))
        target.write('\n')
    target.close()
    print("Bingo", "*" * 50 )
    return GmRs

def Aeq(zeta, Rj, Rk):
    l = 8.00
    n = 4.0
    GmAs = []
    for a in angular_rs_values:
        for b in ts_values:
            GmA = 0
            i = 0
            for e in zeta:
                fk = cutf(Rj[i])
                fj = cutf(Rk[i])
                GmA = GmA + ((( 1+math.cos(e - b)) )** l) * \
                    math.exp(- n *( (((Rj[i] + Rk[i]) / 2) - a) ** 2)) * fj * fk
                i = i+1
            GmA = GmA * (2 ** (1-l))
            GmAs.append(GmA)
    return GmAs

def Apart(atom_types, file_name):
    Fake = Atome_classify()
    C = Fake[0]
    H = Fake[1]
    O = Fake[2]
    target = open(file_name, 'w')
    for a1 in hierarchy.atoms():
        GmRs = []
        e1 = a1.element.upper()
        if e1 in atom_types:
            X = e1
            target.write(str(X))
            target.write('\n')
#fnid angle HH
            Rj = []
            Rk = []
            zeta = []
            for a2 in H:
                if a1 != a2:
                    for a21 in H:
                        if a2 != a21 and a21 != a1:
                            R = a1.distance(a2)
                            Rj.append(R)
                            R = a1.distance(a21)
                            Rk.append(R)
                            ts = a1.angle(a2, a21, deg=False)
                            if ts == None:
                                ts = 0
                            zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('HH')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
#find angle HC
            Rj = []
            Rk = []
            zeta = []
            for a2 in H:
                for a21 in C:
                    if a1 != a2 and a1 != a21 and a2 != a21:
                        R = a1.distance(a2)
                        Rj.append(R)
                        R = a1.distance(a21)
                        Rk.append(R)
                        ts = a1.angle(a2, a21, deg=False)
                        if ts == None:
                            ts = 0
                        zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('HC')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
#find angle HO
            Rj = []
            Rk = []
            zeta = []
            for a2 in H:
                for a21 in O:
                    if a2 != a21 and a21 != a1 and a1 != a2:
                        R = a1.distance(a2)
                        Rj.append(R)
                        R = a1.distance(a21)
                        Rk.append(R)
                        ts = a1.angle(a2, a21, deg=False)
                        if ts == None:
                            ts = 0
                        zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('HO')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
#find angle CC
            Rj = []
            Rk = []
            zeta = []
            for a2 in C:
                for a21 in C:
                    if a2 != a21 and a21 != a1 and a1 != a2:
                        R = a1.distance(a2)
                        Rj.append(R)
                        R = a1.distance(a21)
                        Rk.append(R)
                        ts = a1.angle(a2, a21, deg=False)
                        if ts == None:
                            ts = 0
                        zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('CC')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
#find angle CO
            Rj = []
            Rk = []
            zeta = []
            for a2 in C:
                for a21 in O:
                    if a2 != a21 and a21 != a1 and a1 != a2:
                        R = a1.distance(a2)
                        Rj.append(R)
                        R = a1.distance(a21)
                        Rk.append(R)
                        ts = a1.angle(a2, a21, deg=False)
                        if ts == None:
                            ts = 0
                        zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('CO')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
#find angle OO
            Rj = []
            Rk = []
            zeta = []
            for a2 in O:
                for a21 in O:
                    if a2 != a21 and a21 != a1 and a1 != a2:
                        R = a1.distance(a2)
                        Rj.append(R)
                        R = a1.distance(a21)
                        Rk.append(R)
                        ts = a1.angle(a2, a21, deg=False)
                        if ts == None:
                            ts = 0
                        zeta.append(ts)
            print(zeta, Rj, Rk)
            GmR = Aeq(zeta, Rj, Rk)
            target.write('OO')
            target.write(str(GmR))
            target.write('\n')
            GmRs.append(GmR)
    target.write(str(GmRs))
    target.close()
    print("Bingo", "*" * 50 )
    return GmRs

# pdb_inp = iotbx.pdb.input(file_name="acid060.pdb")
# hierarchy = pdb_inp.construct_hierarchy()
# Rpart(atom_type, file_name='text1.txt')
# Apart(atom_type, file_name='text3.txt')
# pdb_inp = iotbx.pdb.input(file_name="acid300.pdb")
# hierarchy = pdb_inp.construct_hierarchy()
# Rpart(atom_type, file_name='text2.txt')
# Apart(atom_type, file_name='text4.txt')

def get_atom_elements(hierarchy):
    rc = []
    for atom in hierarchy.atoms():
        if atom.element.strip() not in rc:
            rc.append(atom.element.strip())
    return rc


def main(pdb_file_name):
    pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    hierarchy = pdb_inp.construct_hierarchy()
    Rpart(hierarchy, file_name='%s_radial.txt' % pdb_file_name.replace('.pdb',''))
    Apart(atom_type, file_name='%s_angle.txt' % pdb_file_name.replace('.pdb',''))

if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))
