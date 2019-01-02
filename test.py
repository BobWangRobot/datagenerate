import math
import iotbx

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
    i = 0
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
            i = i + 1
    print(GmAs, i)
    return GmAs

i = int(input("inputi"))
j = int(input("inputj"))
a = D2(i, j)
zetas = [1, 2, 3, 4]
l = [0,1]
Rj = zetas
Rs = [1, 2, 3, 4, 5, 6, 7, 8]
Gmr = Rpart(j, 4, Rj, Rs, 6)
#Gma = Apart(i, j, l, a, zetas, 2, 4, 5, 3, 6)
print('g', Gmr)
print(a[1][2])