from iotbx import pdb
import iotbx
import math

class AEV(object):

    def __init__(self, pdb_file_name):
        self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
        self.hierarchy = self.pdb_inp.construct_hierarchy()
        self.radial_nu = 8
        self.radial_angular_nu = 8
        self.eta = 4.0

    def Atome_type(self):
        atom_elements = list()
        for a in self.hierarchy.atoms():
            e = a.element.upper().strip()
            atom_elements.append(e)
        atom_elements = set(atom_elements)
        return atom_elements

    def Atome_classify(self):
        atom_elements = {}
        for a in self.hierarchy.atoms():
            e = a.element.upper().strip()
            atom_elements.setdefault(e, [])
            atom_elements[e].append(a)
        return atom_elements


class Rpart(AEV):
    def __init__(self, pdb_file_name):
        super(Rpart, self).__init__(pdb_file_name)
        self.rs_values = [ 0.900000,1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
        self.Rj =[2.1, 2.2, 2.5]
        self.radial_cutoff = 5.2

    def cutf(R_distance, R_cutoff):
        if R_distance >= R_cutoff:
            Fc = 0.5 * math.cos(math.pi * R_distance / R_cutoff) + 0.5
        else:
            Fc = 0
        return Fc

    def Req(self, Rj):
        GmRs = []
        n = 4.0
        for b in self.rs_values:
            print(self.rs_values)
            GmR = 0
            for a in Rj:
                f = a.cutf(b)
                GmR += math.exp(- n * ((a - b) ** 2)) * f
            GmRs.append(GmR)
        return GmRs

    def Rpart(self):
        assert 0
        atom_types = get_atom_elements(self.hierarchy)
        Fake = Atome_classify(hierarchy)
        for a1 in hierarchy.atoms():
            GmRs = []
            e1 = a1.element.upper().strip()
            assert e1 in atom_types, '%s not in %s' % (e1, atom_types)
            for e2, atoms2 in Fake.items():
                Rj = []
                for a2 in atoms2:
                    R = a1.distance(a2)
                    Rj.append(R)
                GmR = Req(Rj)
                GmRs.append(GmR)
            assert 0
            Rj = []
            for a3 in C:
                R = a1.distance(a3)
                Rj.append(R)
            GmR = Req(Rj)
            GmRs.append(GmR)
            Rj = []
            for a4 in O:
                R = a1.distance(a4)
                Rj.append(R)
            GmR = Req(Rj)
            GmRs.append(GmR)
        print("Bingo", "*" * 50)
        return GmRs

class Apart():
    def __init__(self,ts_values, angular_rs_values, angular_cutoff, angular_zeta):
        super().__init__()
        self.ts_values = ts_values
        self.angular_rs_values = angular_rs_values
        self.angular_cutoff = angular_cutoff
        self.angular_zeta = angular_zeta
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
                    GmA = GmA + (((1 + math.cos(e - b))) ** l) * \
                          math.exp(- n * ((((Rj[i] + Rk[i]) / 2) - a) ** 2)) * fj * fk
                    i = i + 1
                GmA = GmA * (2 ** (1 - l))
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
                # fnid angle HH
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
                # find angle HC
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
                # find angle HO
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
                # find angle CC
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
                # find angle CO
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
                # find angle OO
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
        print("Bingo", "*" * 50)
        return GmRs
    pass

