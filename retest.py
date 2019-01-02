from __future__ import division
import iotbx.pdb
import iotbx.cif
from libtbx import group_args
import mmtbx.model
from libtbx.utils import null_out

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

def read(model):
    Rs = []
    R = 0
    hierarchy = model.get_hierarchy()
    for a1 in hierarchy.atoms():
        for a2 in hierarchy.atoms():
            R = a1.distance(a2)
            Rs.append(R)
            R = 0
            print("kk")
            print(Rs)


def exercise():
    files = [["1qwl.pdb", None]]
    for (pdb_file_name, cif_file_name) in files:
        print(pdb_file_name, "-" * 50)
        model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
        result = read(model=model)

if __name__ == '__main__':
    exercise()