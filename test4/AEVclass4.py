import numpy as np
import iotbx
import math
import time
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
      outl += '  %s :' % (key)
      for v in item:
        outl += '%0.4f, ' % v
      outl += '\n'
    return outl

class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for key1,value in item.items():
        outl += ' %s: '%key1
        outl += '%0.2f, ' % value
      outl += '\n'
    return outl

class AEV(object):

  def __init__(self, pdb_file_name=None, raw_records=None):
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
    self.rs_values = [2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0]
    self.Rj = [2.1, 2.2, 2.5]
    self.cutoff = 8.1
    self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.angular_rs_values = [3.8, 5.2, 5.5, 6.2]
    self.angular_zeta = 8
    self.EAEVs = radial_aev_class()
    self.MAEVs = radial_aev_class()
    self.BAEVs = radial_aev_class()

  # generate ca list
  def generate_ca(self):
    t0 = time.time()
    self.hierarchy.reset_atom_i_seqs()
    self.hierarchy.reset_i_seq_if_necessary()
    for five in generate_protein_fragments(self.hierarchy,
                                           self.geometry_restraints_manager,
                                           include_non_standard_peptides=True,
                                           length=5):
      rc = []
      res_name = []
      for residue in five:
        res_name.append(residue.resname + residue.resseq)
        for atom in residue.atoms():
          if atom.name == ' CA ':
            rc.append(atom)
      if len(rc) == 5:
        rc.append(res_name)
        yield rc
    print('time', time.time() - t0)

  # generate a dictionary about all atome
  # def Atome_classify(self, atype):
  #   atom_elements = {}
  #   for b in self.hierarchy.atoms():
  #     #e = b.element.upper().strip()
  #     e = b.name.strip()
  #     if e == atype:
  #       atom_elements.setdefault(e, [])
  #       atom_elements[e].append(b)
  #   return atom_elements

  def generate_AEV(self):
    for atomlist in self.generate_ca():
      self.center_atom = atomlist[0]
      self.res_name = atomlist[-1][0]
      self.BAEVs.update(self.calculate(atomlist[0:-1]))
      self.center_atom = atomlist[-2]
      self.res_name = atomlist[-1][-1]
      self.EAEVs.update(self.calculate(atomlist[0:-1]))
      self.center_atom = atomlist[2]
      self.res_name = atomlist[-1][2]
      self.MAEVs.update(self.calculate(atomlist[0:-1]))
    return 0

  # cutoff function
  def cutf(self, distance):
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

  def calculate(self, atomlist):
    n = 4.0
    l = 8.0
    AEVs = radial_aev_class()
    res_name = self.res_name
    AEVs.setdefault(res_name, [])
    atom1 = self.center_atom
    atomlist.remove(atom1)
    for Rs in self.rs_values:
      GmR = 0
      for atom2 in atomlist:
        R = atom1.distance(atom2)
        f = self.cutf(R)
        if f != 0:
          mR = math.exp(- n * ((R - Rs) ** 2)) * f
          GmR += mR
      AEVs[res_name].append(GmR)
    for Rs in self.angular_rs_values:
      for zetas in self.ts_values:
        GmA = 0
        for atom2 in atomlist:
          atomlist.remove(atom2)
          for atom3 in atomlist:
            Rij = atom1.distance(atom2)
            Rik = atom1.distance(atom3)
            ZETAijk = atom1.angle(atom2, atom3)
            if ZETAijk != 0:
              fk = self.cutf(Rik)
              fj = self.cutf(Rij)
              if fk != 0 and fj != 0:
                mA = (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                     math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                GmA += mA
        GmA = GmA * 2
        AEVs[res_name].append(GmA)
    return AEVs
