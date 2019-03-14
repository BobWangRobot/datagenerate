import iotbx
from AEVclass import AEV
from AEVclass import Rpart
from AEVclass import Apart

def main(pdb_file_name):
    #a = AEV(pdb_file_name)
    #print(a.Atome_classify())
    #b = Rpart(pdb_file_name)
    #print(b.R_AEV())
    c = Apart(pdb_file_name)
    print(c.A_AEV())

if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

#a = AEV('acid060.pdb')
#print(a.Atome_classify())
#b = Rpart('acid060.pdb')
#print(b.R_AEV())
#c = Apart('acid060.pdb')
#print(c.Aeq())

