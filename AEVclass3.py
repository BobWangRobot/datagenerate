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
        outl += '     %s : %0.3f ' % (e, vec)
      outl += '\n'
    return outl

class AEV(object):

  def __init__(self, direction, scope, pdb_file_name=None, raw_records=None):
    if pdb_file_name:
      self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    else:
      self.pdb_inp = iotbx.pdb.input(lines=raw_records, source_info='raw_records')
    self.direction = direction
    self.hierarchy = self.pdb_inp.construct_hierarchy()
    self.mon_lib_srv = server.server()
    self.ener_lib = server.ener_lib()
    self.processed_pdb = pdb_interpretation.process(self.mon_lib_srv, self.ener_lib, file_name=pdb_file_name,
                                                    raw_records=raw_records)
    self.geometry_restraints_manager = self.processed_pdb.geometry_restraints_manager()
    self.atom_elements = {}
    # self.rs_values = [0.900000, 2.100000, 3.300000, 4.500000, 5.700000, 6.900000, 8.100000, 9.300000]
    self.rs_values = [2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0]
    self.Rj = [2.1, 2.2, 2.5]
    self.cutoff = float(scope)
    # self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.ts_values = [1.178097, 2.748894]
    self.angular_rs_values = [3.8, 5.2, 5.5, 6.2]
    self.angular_zeta = 8
    self.five = []
    self.rc = []
    self.AEVs = radial_aev_class()
    self.FAEVs = radial_aev_class()
    self.BAEVs = radial_aev_class()
    self.diffs = diff_class()
    self.atom_range = diff_class()

  # generate five ca
  def generate_ca(self):
    t0 = time.time()
    self.hierarchy.reset_atom_i_seqs()
    self.hierarchy.reset_i_seq_if_necessary()
    for five in generate_protein_fragments(self.hierarchy,
                                           self.geometry_restraints_manager,
                                           include_non_standard_peptides=True,
                                           length=5):
      rc = []
      for atom in five.atoms():
        if atom.name == ' CA ':
          rc.append(atom)
      if len(rc) == 5:
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
      #e = b.element.upper().strip()
      e = b.name.strip()
      if e == atype:
        atom_elements.setdefault(e, [])
        atom_elements[e].append(b)
    return atom_elements


  def Rpart(self):
    n = 4.0
    AEVs = radial_aev_class()
    atom_range = diff_class()
    dis = self.Atome_classify('CA')
    for atom1 in self.five:
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a + x, {})
      atom_range.setdefault(a + x, {})
      atom_range[a + x].setdefault('Rpart')
      for b, atom2list in dis.items():
        AEVs[a + x].setdefault(b, [])
        for Rs in self.rs_values:
          FGmR = 0
          BGmR = 0
          GmR = 0
          F_range = 0
          B_range = 0
          All_range = 0
          for atom2 in atom2list:
            z = str(atom2.i_seq)
            if atom1 != atom2:
              R = atom1.distance(atom2)
              f = self.cutf(R)
              if f != 0:
                mR = math.exp(- n * ((R - Rs) ** 2)) * f
                GmR += mR
                All_range = All_range + 1
                if int(z) < int(x):
                  FGmR += mR
                  F_range = F_range + 1
                elif int(z) > int(x):
                  BGmR += mR
                  B_range = B_range + 1
          if self.direction=='back':
            AEVs[a+x][b].append(BGmR)
            atom_range[a + x]['Rpart'] = B_range
          if self.direction=='fward':
            AEVs[a+x][b].append(FGmR)
            atom_range[a + x]['Rpart'] = F_range
          if self.direction=='all':
            AEVs[a + x][b].append(GmR)
            atom_range[a + x]['Rpart'] = All_range
    self.AEVs.update(AEVs)
    self.atom_range.update(atom_range)
    return AEVs

  def Apart(self):
    l = 8.00
    n = 4.0
    AEVs = radial_aev_class()
    atom_range = diff_class()
    for atom1 in self.five:
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a+x, {})
      atom_range.setdefault(a + x, {})
      atom_range[a + x].setdefault('Apart')
      for b, atomlist in self.Atome_classify("CA").items():
        for Rs in self.angular_rs_values:
          for zetas in self.ts_values:
            AEVs[a+x].setdefault(b+b, [])
            GmA = 0
            FGmA = 0
            BGmA = 0
            F_repeat_list = []
            B_repeat_list = []
            repeat_list = []
            for atom2 in atomlist:
              y = str(atom2.i_seq)
              for atom3 in atomlist:
                if not atom3.i_seq in repeat_list:
                  z = str(atom3.i_seq)
                  if atom2 != atom1 and atom3 != atom1 and atom2 != atom3:
                    Rij = atom1.distance(atom2)
                    Rik = atom1.distance(atom3)
                    ZETAijk = atom1.angle(atom2, atom3)
                    if ZETAijk !=0:
                      fk = self.cutf(Rik)
                      fj = self.cutf(Rij)
                      if fk != 0 and fj != 0:
                        mA = (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                              math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                        GmA += mA
                        repeat_list.append(y)
                        if int(y) < int(x) and int(z) < int(x):
                          FGmA += mA
                          F_repeat_list.append(y)
                        elif int(y) > int(x) and int(z) > int(x):
                          BGmA += mA
                          B_repeat_list.append(y)
            GmA = GmA * (2 ** (1 - l))
            BGmA = BGmA * (2 ** (1 - l))
            FGmA = FGmA * (2 ** (1 - l))
            if GmA < 1e-6:
              GmA = 0
            if self.direction == 'back':
              AEVs[a + x][b + b].append(BGmA)
              atom_range[a + x]['Apart'] = len(set(B_repeat_list))
            if self.direction == 'fward':
              AEVs[a + x][b + b].append(FGmA)
              atom_range[a + x]['Apart'] = len(set(F_repeat_list))
            if self.direction == 'all':
              AEVs[a + x][b + b].append(GmA)
              atom_range[a + x]['Apart'] = len(set(repeat_list))
      self.AEVs[a + x].update(AEVs[a + x])
      self.atom_range[a + x].update(atom_range[a + x])
    return AEVs

  def get_AEVs(self):
    self.Rpart()
    self.Apart()
    return self.AEVs



