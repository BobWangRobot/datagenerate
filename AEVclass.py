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

#generate a set about all atom types
    def Atome_type(self):
        atom_elements = list()
        for a in self.hierarchy.atoms():
            e = a.element.upper().strip()
            atom_elements.append(e)
        atom_elements = set(atom_elements)
        return atom_elements

#generate a dictionary about all atomes
    def Atome_classify(self):
        atom_elements = {}
        for a in self.hierarchy.atoms():
            e = a.element.upper().strip()
            atom_elements.setdefault(e, [])
            atom_elements[e].append(a)
        return atom_elements

#generate a dictionary include {element_name:distance}
    def All_dis(self):
        All_distance = {}
        for a in self.Atome_type():
            for a1, atom1 in self.Atome_classify().items():
                if a1 == a:
                    for a2, atom2 in self.Atome_classify():
                        if a1 != a2:
                            R = atom1.distace(atom2)
                            All_distance.setdefault(a1, [])
                            All_distance[a1].append(R)
        return All_distance

#generate a dictionary about all angels {element_name: angel}
    def All_angel(self):
        All_angels = {}
        for a in self.Atome_type():
            for a1, atom1 in self.Atome_classify():
                if a1 == a:
                    for b, batom in self.Atome_classify():
                        f = self.Atome_classify()
                        for c, catom in f:
                            angel = atom1.angel(batom, catom)
                            All_angels.setdefault(a1+b+c, [])
                            All_angels[a1+b+c].append(angel)
                        f.pop(b)
        return All_angels

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





class Apart(AEV):
    def __init__(self, pdb_file_name):
        super(Rpart, self).__init__(pdb_file_name)
        self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
        self.angular_rs_values = [0.900000, 2.200000]
        self.angular_cutoff = 3.5
        self.angular_zeta = 8

    def Aeq(self, Rj, Rk):
        l = 8.00
        n = 4.0
        GmAs = []
        for a in self.angular_rs_values:
            for b in self.ts_values:
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

    def Apart(self):
        Fake = Atome_classify()
        C = Fake[0]
        H = Fake[1]
        O = Fake[2]
        for a1 in hierarchy.atoms():
            GmRs = []
            e1 = a1.element.upper()
            if e1 in atom_types:
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
                GmRs.append(GmR)
        print("Bingo", "*" * 50)
        return GmRs
    pass

