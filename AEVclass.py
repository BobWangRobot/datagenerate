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
                    All_distance.setdefault(a1, [])
                    for b in self.Atome_type():
                        for b1, batom1 in self.Atome_classify().items():
                            if b1 == b:
                                All_distance[a1].setdefault(b1, [])
                                R = atom1.distace(atom2)
                                All_distance[a1:b1].append(R)
        return All_distance

#generate a dictionary about all angels {element_name: angel}
    def All_angel(self):
        All_angels = {}
        for a in self.Atome_type():
            for a1, atom1 in self.Atome_classify().items():
                if a1 == a:
                    for b, batom in self.Atome_classify().items():
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

    def cutf(self,R_distance):
        if R_distance >= self.radial_cutoff:
            Fc = 0.5 * math.cos(math.pi * R_distance / self.radial_cutoff) + 0.5
        else:
            Fc = 0
        return Fc

    def R_AEV(self, Rj):
        AEVs = {}
        n = 4.0
        for a, atom1 in self.All_dis().items():
            AEVs.setdefaut(a, [])
            for Rs in self.rs_values:
                GmR = 0
                for R in atom1:
                    f = R.cutf(self.radial_cutoff)
                    GmR += math.exp(- n * ((R - Rs) ** 2)) * f
                AEVs[a].append(GmR)
        return AEVs

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


    pass

