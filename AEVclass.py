from iotbx import pdb
import numpy
import iotbx
import math

class AEV(object):

  def __init__(self, pdb_file_name):
    self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    self.hierarchy = self.pdb_inp.construct_hierarchy()
    self.radial_nu = 8
    self.radial_angular_nu = 8
    self.eta = 4.0
#cutoff function
  def cutf(self, distance):
      if distance <= self.cutoff:
        Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
      else:
        Fc = 0
      return Fc

#generate a dictionary about all atomes
  def Atome_classify(self):
      atom_elements = {}
      list1 = list()
      for a in self.hierarchy.atoms():
          e = a.element.upper().strip()
          list1.append(e)
      list1 = set(list1)
      for b in self.hierarchy.atoms():
        e = b.element.upper().strip()
        atom_elements.setdefault(e, [])
        atom_elements[e].append(b)
      return atom_elements




class Rpart(AEV):
    def __init__(self, pdb_file_name):
        super(Rpart, self).__init__(pdb_file_name)
        self.rs_values = [0.900000, 1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
        self.Rj =[2.1, 2.2, 2.5]
        self.cutoff = 3.5


    def R_AEV(self):
        AEVs = {}#It is a multiple dictionary
        n = 4.0
        for a, atom_list1 in self.Atome_classify().items():
            AEVs.setdefault(a, {})
            for Rs in self.rs_values:
                for atom1 in atom_list1:
                    for b, atom_list2 in self.Atome_classify().items():
                        AEVs[a].setdefault(b, [])
                        GmR = 0
                        for atom2 in atom_list2:
                            R = atom1.distance(atom2)
                            f = self.cutf(R)
                            GmR += math.exp(- n * ((R - Rs) ** 2)) * f
                        AEVs[a][b].append(GmR)
        return AEVs

class Apart(AEV):
    def __init__(self, pdb_file_name):
        super(Rpart, self).__init__(pdb_file_name)
        self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
        self.angular_rs_values = [0.900000, 2.200000]
        self.angular_cutoff = 3.5
        self.angular_zeta = 8
        self.radial_cutoff = 5.2

    def Aeq(self):
        l = 8.00
        n = 4.0
        GmAs = []
        for a in self.angular_rs_values:
            for b in self.ts_values:
                GmA = 0
                i = 0
                for e in zeta:
                    fk = self.angular_cutoff.cutf()
                    fj = self.radial_cutoff.cutf()
                    GmA = GmA + (((1 + math.cos(e - b))) ** l) * \
                          math.exp(- n * ((((Rj[i] + Rk[i]) / 2) - a) ** 2)) * fj * fk
                    i = i + 1
                GmA = GmA * (2 ** (1 - l))
                GmAs.append(GmA)
        return GmAs
