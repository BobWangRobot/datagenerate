import numpy as np
import iotbx
import math
import time
import copy
from itertools import chain
from collections import OrderedDict
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
          outl += '%0.3f, ' % v
        outl += '\n'
    return outl

class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :\n' % (key)
      for e, vec in item.items():
        outl += '     %s : %0.2f ' % (e, vec)
      outl += '\n'
    return outl

class AEV_base(object):

  def __init__(self, pdb_file_name=None, raw_records=None):
    # assert count(pdb_file_name, raw_records)==1
    if pdb_file_name:
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
    self.rc = []

  # generate five ca
  def generate_ca(self):
    t0 = time.time()
    self.hierarchy.reset_atom_i_seqs()
    self.hierarchy.reset_i_seq_if_necessary()
    for five in generate_protein_fragments(self.hierarchy,
                                           self.geometry_restraints_manager,
                                           include_non_standard_peptides=True,
                                           length=2):
      rc = []
      for atom in five.atoms():
        if atom.name == ' CA ':
          rc.append(atom)
      if len(rc) == 2:
        yield rc
    print 'time', time.time() - t0

  # cutoff function
  def cutf(self, distance):
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

  # generate a dictionary about all atome
  def Atome_classify(self, atype):
    atom_elements = {}
    for b in self.hierarchy.atoms():
      e = b.element.upper().strip()
      #e = b.name.strip()
      if e == atype:
        atom_elements.setdefault(e, [])
        atom_elements[e].append(b)
    return atom_elements


class AEV(AEV_base):
  def __init__(self, pdb_file_name=None, raw_records=None):
    AEV_base.__init__(self, pdb_file_name, raw_records)
    self.rs_values = [0.900000, 2.100000, 3.300000, 4.500000, 5.700000, 6.900000, 8.100000, 9.300000]
    self.Rj = [2.1, 2.2, 2.5]
    self.radial_cutoff = 10.5
    self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.angular_rs_values = [0.900000, 6.225000]
    self.angular_cutoff = 8
    self.angular_zeta = 8
    self.five = []
    self.AEVs = radial_aev_class()
    self.rdistance = radial_aev_class()
    self.diffs = diff_class()

  def Rpart(self):
    n = 4.0
    AEVs = radial_aev_class()
    dis = self.Atome_classify('C')
    for atom1 in self.five:
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a + x, {})
      for b, atom2list in dis.items():
        AEVs[a + x].setdefault(b, [])
        for Rs in self.rs_values:
          GmR = 0
          for atom2 in atom2list:
            if atom1 != atom2:
              R = atom1.distance(atom2)
              self.cutoff = self.radial_cutoff
              f = self.cutf(R)
              if f != 0:
                GmR += math.exp(- n * ((R - Rs) ** 2)) * f
          AEVs[a+x][b].append(GmR)
          self.AEVs.update(AEVs)
    return AEVs

  def Apart(self):
    l = 8.00
    n = 4.0
    AEVs = radial_aev_class()
    for atom1 in self.five:
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a+x, {})
      f = dict(self.Atome_classify("C"))
      for b, atom2list in self.Atome_classify("C").items():
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
        self.AEVs[a+x].update(AEVs[a+x])
    return AEVs

  def get_AEVs(self):
    self.Rpart()
    self.Apart()
    return self.AEVs

  def get_items(self):
    empty = self.get_AEVS()
    for ele, values in empty.items():
      for item in values.keys():
        empty[ele][item] = [0, 0, 0, 0, 0, 0, 0, 0]
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


  def compare(self, match):
    for ele1,item1 in self.Rpart().items():
      self.diffs.setdefault(ele1, OrderedDict())
      for ele2,item2 in match.Rpart().items():
        for list1,list2 in zip(item1.values(),item2.values()):
          covalue = np.corrcoef(list1, list2).tolist()
          if covalue[1][0] < 0.9:
            covalue[1][0] = 0
          self.diffs[ele1].setdefault(ele2, covalue[1][0])
      # all += covalue[1][0]
      # if covalue[1][0] > limit:
      #   self.diffs.setdefault(ele1 + ele2, covalue[1][0])

  def find_function(self, match):
    for self.five in self.generate_ca():
      for match.five in match.generate_ca():
      #for match.five in match.generate_ca():
        self.compare(match)
    print(self.diffs)
        #print(diffs)
        # if diffs['all'] > 0.9:
        #   print(diffs.keys())

