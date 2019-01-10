from __future__ import division
import math
import iotbx.pdb
import iotbx.cif
from libtbx import group_args
import mmtbx.model
from libtbx.utils import null_out

# AEVs value
rs_values = [0.900000, 1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
angular_rs_values = [0.900000, 2.200000]
radial_cutoff = 5.2
angular_cutoff = 3.5
radial_nu = 8
radial_angular_nu = 8
angular_zeta = 8

##get model of pdb file
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

##read model and get information of distance and angle
def read(model, ele_name):
    Rs = []
    zeta = []
    hierarchy = model.get_hierarchy()
    for e in ele_name:
        for a1 in hierarchy.atoms():
            e1 = a1.element.upper()
            if e1 == e:
                for a2 in hierarchy.atoms():
                    e2 = a2.element.upper
                    while a2 != a1:
                        Rj = a1.distance(a2)
                        Ris.append(R)
                        for a3 in hierarchy.atoms():
                            while a3 != a2:
                                A = a1.angle(a2, a3, deg = True)
                                if A == None:continue
                                zeta.append(A)

    return e, Rs, zeta


def cutf(Rij,Rc):
    if Rc >= Rij:
       Fc = 0.5*math.cos(math.pi*Rij/Rc)+0.5
    else:
       Fc = 0
    return Fc

def Rpart(Rc, Rj, Rs, n):
    GmRs = []
    for b in Rs:
        GmR = 0
        for a in Rj:
            f = cutf(a, Rc)
            GmR = GmR + math.exp(n * ((a - b) ** 2))*f
            GmRs.append(GmR)
    return GmRs

def Apart(l, zeta, zetas, n, Rj, Rk, Rs, Rc):
    l = 8.00
    GmAs = []
    fj = cutf(Rj, Rc)
    fk = cutf(Rk, Rc)
    for d in Rs:
        for c in zetas:
            GmA = 0
            for a in zeta:
                GmA = GmA + (((1+math.cos(a - c)) )** l) * \
                          math.exp(-n *( ((Rj + Rk) / 2 - d) ** 2)) * fj * fk
            GmA = GmA + 2 ** (1-l)
            GmAs.append(GmA)
    return GmAs


def exercise():
    files = [["1qwm.pdb", None]]
    for (pdb_file_name, cif_file_name) in files:
        print(pdb_file_name, "-" * 50)
        model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
        Inf = read(model=model)
        Element_namr = Inf[0]
        Dis = Inf[1]
        Ang = Inf[2]
        Rnum = Rpart(radial_cutoff, Dis, rs_values, 4.0)
        Rv = Rpart(radial_cutoff, Rs, rs_values, 4.0)
        Av = Apart(l, zeta, ts_values, 4.0, Rjs, Rks, rs_values)

        Anum =
        print(result)


if __name__ == '__main__':
    exercise()
