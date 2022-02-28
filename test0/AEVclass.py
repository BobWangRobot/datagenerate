from iotbx import pdb
import numpy as np
import iotbx
import math
from collections import OrderedDict

class radial_aev_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :\n' % (key)
      for e, vec in item.items():
        outl += '     %s : ' % e
        for v in vec:
          outl += '%0.3f,' % v
        outl += '\n'
    return outl


class AEV_base(object):

  def __init__(self, pdb_file_name=None, raw_records=None):
    #assert count(pdb_file_name, raw_records)==1
    if pdb_file_name :
      self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    else:
      self.pdb_inp = iotbx.pdb.input(lines=raw_records, source_info='raw_records')
    self.hierarchy = self.pdb_inp.construct_hierarchy()

#cutoff function
  def cutf(self, distance):
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

#generate a dictionary about all atome
  def Atome_classify(self):
    atom_elements = {}
    for b in self.hierarchy.atoms():
      e = b.element.upper().strip()
      atom_elements.setdefault(e, [])
      atom_elements[e].append(b)
    return atom_elements

class AEV(AEV_base):
  def __init__(self, pdb_file_name=None, raw_records=None):
    AEV_base.__init__(self, pdb_file_name, raw_records)
    self.rs_values = [0.900000, 1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
    self.Rj = [2.1, 2.2, 2.5]
    self.radial_cutoff = 5.2
    self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.angular_rs_values = [0.900000, 2.200000]
    self.angular_cutoff = 3.5
    self.angular_zeta = 8

  def Rpart(self):
    AEVs = radial_aev_class() #It is a multiple dictionary
    n = 4.0
    for atom1 in self.hierarchy.atoms():
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a+x, OrderedDict())
      for b, atom2list in self.Atome_classify().items():
        for Rs in self.rs_values:
          AEVs[a+x].setdefault(b, [])
          GmR = 0
          for atom2 in atom2list:
            if atom1 != atom2:
              R = atom1.distance(atom2)
              self.cutoff = self.radial_cutoff
              f = self.cutf(R)
              if f != 0:
                GmR += math.exp(- n * ((R - Rs) ** 2)) * f
              else: continue
          if GmR<1e-6:
            GmR = 0
          AEVs[a+x][b].append(GmR)
    return AEVs

  def Apart(self):
    l = 8.00
    n = 4.0
    AEVs = radial_aev_class()
    for atom1 in self.hierarchy.atoms():
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a+x, OrderedDict())
      f = self.Atome_classify()
      for b, atom2list in f.items():
        for c, atom3list in f.items():
          for Rs in self.angular_rs_values:
            for zetas in self.ts_values:
              AEVs[a+x].setdefault(b+c, [])
              GmA = 0
              for atom2 in atom2list:
                for atom3 in atom3list:
                  if atom2 != atom1 and atom3 != atom1:
                    Rij = atom1.distance(atom2)
                    Rik = atom1.distance(atom2)
                    ZETAijk = atom1.angle(atom2, atom3)
                    self.cutoff = self.angular_cutoff
                    fk = self.cutf(Rik)
                    fj = self.cutf(Rij)
                    if fk != 0 and fj != 0:
                      GmA += (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                            math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                    else: continue
              GmA = GmA * (2 ** (1 - l))
              if GmA < 1e-6:
                GmA = 0
              AEVs[a+x][b+c].append(GmA)
        f.pop(b)#delecte repeated atomes
    return AEVs

  def get_AEVS(self):
    all_AEV = self.Rpart()
    for element, Rvc in all_AEV.items():
      Rvc.update(self.Apart()[element])
    return all_AEV

  def get_items(self):
    empty = self.get_AEVS()
    for ele, values in empty.items():
      for item in values.keys():
        empty[ele][item] = [0,0,0,0,0,0,0,0]
    return empty

  def merge(self, b):
    a = self.get_AEVS()
    for key, item in b.items():
      if key in a:
        for key1, item1 in item.items():
          if key1 not in a[key].keys():
            a[key].setdefault(key1, [])
            a[key][key1] = item1
      else:
        a.setdefault(key, {})
        a[key] = item
    return a

  def compare(self, match_item, element_list=None):
    aev1 = self.merge(match_item.get_items())
    aev2 = match_item.merge(self.get_items())
    print(aev1, aev2)
    diff = {}
    if element_list:
      list = element_list
    else:
      list = aev2.keys()
    for element in list:
      diff.setdefault(element, [])
      all1 = []
      for r_or_a, value in aev1[element].items():
        for v1 in value:
          all1.append(v1)
      for element2 in list:
        all2 = []
        try:
          for r_or_a2 in aev1[element].keys():
            value2 = aev2[element2][r_or_a2]
            for v2 in value2:
              all2.append(v2)
          covalue = np.corrcoef(all1, all2).tolist()
          diff[element].append(covalue[1][0])
        except KeyError:
          print('%s:type error'%element)
          continue
    return diff

  # def copare_detail(self,match):
