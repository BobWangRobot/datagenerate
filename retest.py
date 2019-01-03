from __future__ import division
import iotbx.pdb
import iotbx.cif
from libtbx import group_args
import mmtbx.model
from libtbx.utils import null_out

ty = ['C', 'H', 'O']
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

def read1(model):
    Rs = []
    R = 0
    hierarchy = model.get_hierarchy()
    for a1 in hierarchy.atoms():
        for a2 in hierarchy.atoms():
            R = a1.distance(a2)
            Rs.append(R)
            R = 0
            print(Rs)

def read(model, ele_name):
    Rjs = []
    Rks = []
    zeta = []
    hierarchy = model.get_hierarchy()
    for e in ele_name:
        for a1 in hierarchy.atoms():
            e1 = a1.element.upper()
            while e1 == e:
                print(e)
                for a2 in hierarchy.atoms():
                    e2 = a2.element.upper
                    while a2 != a1:
                        Rj = a1.distance(a2)
                        Rjs.append(Rj)
                        for a3 in hierarchy.atoms():
                            while a3 != a2:
                                Rk = a2.distance(a3)
                                A = a1.angle(a2, a3, deg = True)
                                Rks.append(Rk)
                                zeta.append(A)
    print(e, Rjs, Rks,zeta)
    return e, Rjs, Rks, zeta



def exercise():
    files = [["1qwl.pdb", None]]
    for (pdb_file_name, cif_file_name) in files:
        print(pdb_file_name, "-" * 50)
        model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
        result = read(model, ty)

if __name__ == '__main__':
    exercise()