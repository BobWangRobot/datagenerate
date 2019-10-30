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

  def __init__(self, scope, pdb_file_name=None, raw_records=None):
    if pdb_file_name:
      self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    else:
      self.pdb_inp = iotbx.pdb.input(lines=raw_records, source_info='raw_records')
    self.scope = scope
    self.hierarchy = self.pdb_inp.construct_hierarchy()
    self.mon_lib_srv = server.server()
    self.ener_lib = server.ener_lib()
    self.processed_pdb = pdb_interpretation.process(self.mon_lib_srv, self.ener_lib, file_name=pdb_file_name,
                                                    raw_records=raw_records)
    self.geometry_restraints_manager = self.processed_pdb.geometry_restraints_manager()
    self.atom_elements = {}
    self.rs_values = [0.900000, 2.100000, 3.300000, 4.500000, 5.700000, 6.900000, 8.100000, 9.300000]
    self.Rj = [2.1, 2.2, 2.5]
    self.radial_cutoff = 10.5
    self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.angular_rs_values = [0.900000, 6.225000]
    self.angular_cutoff = 8
    self.angular_zeta = 8
    self.five = []
    self.rc = []
    self.AEVs = radial_aev_class()
    self.FAEVs = radial_aev_class()
    self.BAEVs = radial_aev_class()
    self.rdistance = radial_aev_class()
    self.diffs = diff_class()

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
    dis = self.Atome_classify('CA')
    for atom1 in self.five:
      x = str(atom1.i_seq)
      a = atom1.element.upper().strip()
      AEVs.setdefault(a + x, {})
      for b, atom2list in dis.items():
        AEVs[a + x].setdefault(b, [])
        for Rs in self.rs_values:
          FGmR = 0
          BGmR = 0
          GmR = 0
          for atom2 in atom2list:
            z = str(atom2.i_seq)
            if atom1 != atom2:
              R = atom1.distance(atom2)
              self.cutoff = self.radial_cutoff
              f = self.cutf(R)
              if f != 0:
                mR = math.exp(- n * ((R - Rs) ** 2)) * f
                GmR += mR
                if int(z) < int(x):
                  FGmR += mR
                elif int(z) > int(x):
                  BGmR += mR
          if self.scope=='back':
            AEVs[a+x][b].append(BGmR)
          if self.scope=='fward':
            AEVs[a+x][b].append(FGmR)
          if self.scope=='all':
            AEVs[a + x][b].append(GmR)
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
      f = dict(self.Atome_classify("CA"))
      for b, atom2list in self.Atome_classify("CA").items():
        for c, atom3list in f.items():
          for Rs in self.angular_rs_values:
            for zetas in self.ts_values:
              AEVs[a+x].setdefault(b+c, [])
              GmA = 0
              FGmA = 0
              BGmA = 0
              for atom2 in atom2list:
                y = str(atom2.i_seq)
                for atom3 in atom3list:
                  z = str(atom3.i_seq)
                  if atom2 != atom1 and atom3 != atom1:
                    Rij = atom1.distance(atom2)
                    Rik = atom1.distance(atom2)
                    ZETAijk = atom1.angle(atom2, atom3)
                    self.cutoff = self.angular_cutoff
                    fk = self.cutf(Rik)
                    fj = self.cutf(Rij)
                    if fk != 0 and fj != 0:
                      mA = (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                            math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                      GmA += mA
                      if int(y) < int(x) and int(z) < int(x):
                        FGmA += mA
                      elif int(y) > int(x) and int(z) > int(x):
                        BGmA += mA
              GmA = GmA * (2 ** (1 - l))
              BGmA = BGmA * (2 ** (1 - l))
              FGmA = FGmA * (2 ** (1 - l))
              if GmA < 1e-6:
                GmA = 0
              if self.scope == 'back':
                AEVs[a + x][b + c].append(BGmA)
              if self.scope == 'fward':
                AEVs[a + x][b + c].append(FGmA)
              if self.scope == 'all':
                AEVs[a + x][b + c].append(GmA)
        f.pop(b)#delecte repeated atomes
      self.AEVs[a+x].update(AEVs[a+x])
    return AEVs

  def get_AEVs(self):
    self.Rpart()
    self.Apart()
    return self.AEVs



