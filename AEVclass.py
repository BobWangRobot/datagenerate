from iotbx import pdb
import numpy
import iotbx
import math

class radial_aev_class(dict):
  def __repr__(self):
    outl = 'AEV ...\n'
    for key, item in self.items():
      outl += '  %s :\n' % (key)
      for e, vec in item.items():
        outl += '     %s : ' % e
        for v in vec:
          outl += '%0.3f,' % v
        outl += '\n'
    return outl

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
      for b in self.hierarchy.atoms():
        e = b.element.upper().strip()
        atom_elements.setdefault(e, [])
        atom_elements[e].append(b)
      return atom_elements



class Rpart(AEV):
    def __init__(self, pdb_file_name):
        self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
        self.hierarchy = self.pdb_inp.construct_hierarchy()
        self.rs_values = [0.900000, 1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
        self.Rj =[2.1, 2.2, 2.5]
        self.cutoff = 3.5

    def R_AEV(self):
        AEVs = radial_aev_class() 
        n = 4.0
        for atom1 in self.hierarchy.atoms():
            x = str(atom1.i_seq)
            a = atom1.element.upper().strip()
            AEVs.setdefault(a+x, {})
            for b, atom2list in self.Atome_classify().items():
                for Rs in self.rs_values:
                    AEVs[a+x].setdefault(b, [])
                    GmR = 0
                    for atom2 in atom2list:
                        if atom1 != atom2:
                            R = atom1.distance(atom2)
                            f = self.cutf(R)
                            GmR += math.exp(- n * ((R - Rs) ** 2)) * f
                    AEVs[a+x][b].append(GmR)
        return AEVs


class Apart(AEV):
    def __init__(self, pdb_file_name):
        self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
        self.hierarchy = self.pdb_inp.construct_hierarchy()
        self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
        self.angular_rs_values = [0.900000, 2.200000]
        self.cutoff = 3.5
        self.angular_zeta = 8

    def A_AEV(self):
        l = 8.00
        n = 4.0
        AEVs = radial_aev_class()
        for atom1 in self.hierarchy.atoms():
            x = str(atom1.i_seq)
            a = atom1.element.upper().strip()
            AEVs.setdefault(a+x, {})
            for b, atom2list in self.Atome_classify().items():
                for c, atom3list in self.Atome_classify().items():
                    for Rs in self.angular_rs_values:
                        for zetas in self.ts_values:
                            AEVs[a+x].setdefault(b+c, [])
                            GmA = 0
                            for atom2 in atom2list:
                                for atom3 in atom3list:
                                    if atom2 != atom1 and atom3 != atom1 and atom2 != atom3:
                                        Rij = atom1.distance(atom2)
                                        Rik = atom1.distance(atom2)
                                        ZETAijk = atom1.angle(atom2, atom3)
                                        fk = self.cutf(Rik)
                                        fj = self.cutf(Rij)
                                        GmA += (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                                              math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                            GmA = GmA * (2 ** (1 - l))
                            AEVs[a+x][b+c].append(GmA)
        return AEVs
