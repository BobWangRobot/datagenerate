from __future__ import division
import iotbx.pdb
import iotbx.cif
from libtbx import group_args
import mmtbx.model
from libtbx.utils import null_out



def is_bonded(atom_1, atom_2, bond_proxies_simple):
    result = False
    for proxy in bond_proxies_simple:
        i_seq, j_seq = proxy.i_seqs
        if (atom_1.i_seq in proxy.i_seqs and atom_2.i_seq in proxy.i_seqs):
            result = True
            break
    return result


def get_model(pdb_file_name, cif_file_name):
    pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    restraint_objects = None
    if (cif_file_name is not None):
        cif_object = iotbx.cif.reader(cif_file_name).model()
        restraint_objects = [(cif_file_name, cif_object)]
    model = mmtbx.model.manager(
        model_input=pdb_inp,
        restraint_objects=restraint_objects,
        process_input=True,
        log=null_out())
    return model


def find_FA(model):
    i = 0
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=model.get_sites_cart())
    hierarchy = model.get_hierarchy()
    for a1 in hierarchy.atoms():
        e1 = a1.element.upper()
        if (e1 == "O"):
            for a2 in hierarchy.atoms():
                e2 = a2.element.upper()
                if (e2 == "C"):
                    if (is_bonded(a1, a2, bond_proxies_simple)):
                        for a3 in hierarchy.atoms():
                            e3 = a3.element.upper()
                            if (e3 == "O"):
                                if (is_bonded(a1, a3, bond_proxies_simple)):
                                    for a4 in hierarchy.atoms():
                                        e4 = a4.element.upper()
                                        if (e4 == "H"):
                                            if (is_bonded(a1, a4, bond_proxies_simple)):
                                                for a5 in hierarchy.atoms():
                                                    e5 = a5.element.upper()
                                                    if (e5 == "H"):
                                                        if (is_bonded(a2, a5, bond_proxies_simple)):
                                                            print("find formic acid")
        else:
            i = i + 1
            print("None", i)


def exercise():
    files = [["1qwm.pdb", None]]
    for (pdb_file_name, cif_file_name) in files:
        print(pdb_file_name, "-" * 50)
        model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
        result = find_FA(model=model)
        

if __name__ == '__main__':
    exercise()
