import numpy as np
import iotbx
import math
import time
import copy
from collections import OrderedDict
from iotbx import pdb
from iotbx import pdb
from mmtbx import *
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.conformation_dependent_library import generate_protein_fragments

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
    self.mon_lib_srv = server.server()
    self.ener_lib = server.ener_lib()
    self.processed_pdb = pdb_interpretation.process(self.mon_lib_srv, self.ener_lib, file_name=pdb_file_name,
                                                    raw_records=raw_records)
    self.geometry_restraints_manager = self.processed_pdb.geometry_restraints_manager()
    self.atom_elements = {}
    
  #generate five ca
  def generate_ca(self):
    t0 = time.time()
    self.hierarchy.reset_atom_i_seqs()
    self.hierarchy.reset_i_seq_if_necessary()
    for five in generate_protein_fragments(self.hierarchy, self.geometry_restraints_manager, length=5):
      rc = []
      for atom in five.atoms():
        if atom.name == ' CA ':
          rc.append(atom)
      if len(rc) == 5:
        yield rc
    print 'time', time.time() - t0


  #cutoff function
  def cutf(self, distance):
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

#generate a dictionary about all atome
  def Atome_classify(self):
    for b in self.five:#.hierarchy.atoms():
      e = b.element.upper().strip()
      self.atom_elements.setdefault(e, [])
      self.atom_elements[e].append(b)
    return self.atom_elements

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
    self.AEVs = radial_aev_class()
    self.five = []

  def get_AEVS(self):
    n = 4.0
    l = 8.00
    i = 0
    for atom1 in self.five:
      i = i + 1
      x = str(i)
      a = atom1.element.upper().strip()
      self.AEVs.setdefault(a+x, {})
      dis = self.Atome_classify()
      for b, atom2list in dis.items():
        for Rs in self.rs_values:
          self.AEVs[a+x].setdefault(b, [])
          GmR = 0
          for atom2 in atom2list:#radial caculation
            if atom1 != atom2:
              R = atom1.distance(atom2)
              self.cutoff = self.radial_cutoff
              f = self.cutf(R)
              if f != 0:
                GmR += math.exp(- n * ((R - Rs) ** 2)) * f
              else: continue
          if GmR<1e-6:
            GmR = 0
          self.AEVs[a+x][b].append(GmR)
        for c, atom3list in dis.items():#angle caculation
          for Rs in self.angular_rs_values:
            for zetas in self.ts_values:
              self.AEVs[a + x].setdefault(b + c, [])
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
                    else:
                      continue
              GmA = GmA * (2 ** (1 - l))
              if GmA < 1e-6:
                GmA = 0
              self.AEVs[a + x][b + c].append(GmA)
        dis.pop(b)  # delecte repeated atomes
    return self.AEVs

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
  
  def compare(self,match):
    aev1 = self.get_AEVS()
    aev2 = match.get_AEVS()
    diff = {}
    print(aev1, aev2)
    for ele,value in aev1.items():
      diff.setdefault(ele, [])
      all = []
      all1 = []
      for a in value.values():
        all.extend(a)
      for b in aev2[ele].values():
        all1.extend(b)
      covalue = np.corrcoef(all, all1).tolist()
      diff[ele].append(covalue[1][0])
    return diff
    
  
  # def compare(self, match_item, element_list=None):
  #   aev1 = self.get_AEVS()
  #   aev2 = match_item.get_items()
  #   diff = {}
  #   if element_list:
  #     list = element_list
  #   else:
  #     list = aev2.keys()
  #   for element in list:
  #     diff.setdefault(element, [])
  #     all1 = []
  #     for r_or_a, value in aev1[element].items():
  #       for v1 in value:
  #         all1.append(v1)
  #     for element2 in list:
  #       all2 = []
  #       try:
  #         for r_or_a2 in aev1[element].keys():
  #           value2 = aev2[element2][r_or_a2]
  #           for v2 in value2:
  #             all2.append(v2)
  #         covalue = np.corrcoef(all1, all2).tolist()
  #         diff[element].append(covalue[1][0])
  #       except KeyError:
  #         print('%s:type error'%element)
  #         continue
  #   return diff

  def find_function(self):
    for self.five in self.generate_ca():
      print(self.five)
      print(self.get_AEVS())
      self.AEVs = radial_aev_class()
      self.five = []
      
