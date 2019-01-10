import math

rs_values = [ 0.900000,1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
angular_rs_values = [0.900000, 2.200000]
atom_type = [' C', ' H',' O']
radial_cutoff = 5.2
angular_cutoff = 3.5
radial_nu = 8
radial_angular_nu = 8
angular_zeta = 8

def read():
    return 0

def D3(i,j,k):
    a = []
    d = 1
    for b in range(0,i):
        tmp = []
        for c in range(0,j):
            app = []
            for e in range(1,k):
                app.append(d)
                d = d + 1
            tmp.append(app)
        a.append(tmp)
    print(a)
    return a

def D2(i,j):
    a = []
    d = 1
    for b in range(0,i):
        tmp = []
        for b in range(0,j):
            d = d + 1
            tmp.append(d)
        a.append(tmp)
    print(a)
    return a

def cutf(R):
    if angular_cutoff >= R:
       Fc = 0.5*math.cos(math.pi*R/radial_cutoff)+0.5
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



def Apart(zeta, Rj, Rk):
    l = 8.00
    zetas = ts_values
    n = 4.0
    Rs = angular_rs_values
    GmAs = []
    for a in Rs:
        for b in zetas:
            GmA = 0
            for c in Rj:
                for d in Rk:
                    if c > 0:
                        if d > 0:
                            fj = cutf(c)
                            fk = cutf(d)
                            for e in zeta:
                                GmA = GmA + (((1+math.cos(e - b)) )** l) * \
                                          math.exp(- n *( (((d + c) / 2) - a) ** 2)) * fj * fk
            GmA = GmA + 2 ** (1-l)
            GmAs.append(GmA)
    return GmAs

a = [1, 2, 3, 4]
b = [0,1]
c = [1, 2, 3, 4]
d = Apart(a, b, c)
print(d)
