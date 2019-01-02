import math
import iotbx.pdb

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
    return Rs

def find_FA(model):
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=model.get_sites_cart())
    hierarchy = model.get_hierarchy()
    for a1 in hierarchy.atoms():

        e1 = a1.element.upper()
        if (e1 == "C"):
            for a2 in hierarchy.atoms():
                e2 = a2.element.upper()
                if (e2 == "H"):
                    if (is_bonded(a1, a2, bond_proxies_simple)):
                        for a3 in hierarchy.atoms():
                            e3 = a3.element.upper()
                            if (e3 == "O"):
                                if (is_bonded(a1, a3, bond_proxies_simple)):
                                    for a4 in hierarchy.atoms():
                                        e4 = a4.element.upper()
                                        if (e4 == "O"):
                                            if (is_bonded(a1, a4, bond_proxies_simple)):
                                                for a5 in hierarchy.atoms():
                                                    e5 = a5.element.upper()
                                                    if (e5 == "H"):
                                                        if (is_bonded(a4, a5, bond_proxies_simple)):
                                                            print("find formic acid")

def cutf(Rij,Rc):
    if Rc >= Rij:
       Fc = 0.5*math.cos(math.pi*Rij/Rc)+0.5
    else:
       Fc = 0
    return Fc

def Rpart(j, Rc, Rj, Rs, n):
    GmRs = []
    for b in range(8):
        GmR = 0
        for a in range(0,j-1):
            f = cutf(Rj[a], Rc)
            GmR = GmR + math.exp(n * ((Rj[a] - Rs[b]) ** 2))*f
            GmRs.append(GmR)
    return GmRs

def Apart(j, k, l, zeta, zetas, n, Rj, Rk, Rs, Rc):
    GmAs = []
    fj = cutf(Rj, Rc)
    fk = cutf(Rk, Rc)
    for d in range(2):
        for c in range(4):
            GmA = 0
            for a in range(0,j-1):
                for b in range(0, k-1):
                    GmA = GmA + (((1+math.cos(zeta[a][b] - zetas[c])) )** l[d]) * \
                          math.exp(-n *( ((Rj + Rk) / 2 - Rs) ** 2)) * fj * fk
            GmA = GmA + 2 ** (1-l[d])
            GmAs.append(GmA)
    return GmAs