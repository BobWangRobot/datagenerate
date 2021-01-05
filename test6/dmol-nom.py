import math
def read_c_f():
  with open('dmol.outmol', 'r') as f:
    coordinate = []
    force = []
    atom = []
    energy = []
    x = []
    y = []
    z = []
    d = []
    i = 0
    for data in f:
      if not data:
        return None
      if data.startswith('df   Si'):
        l = data.split()
        atom.append(l[1])		
      if data.startswith('df   Si'):
        l = data.split()
        coordinate.append(l[2:5])
        force.append(l[5:])
      elif data.startswith('opt=='):
        l = data.split()
        if len(l) == 6:
          e = l[2]
          energy.append(float(l[2]))
          print (e)
#-------find minimum distance-----------------------------------------------------------------
#        for k in range(len(atom)):
#          for l in range(len(atom)):
#            x.append(math.pow(float(coordinate[k][0])-float(coordinate[l][0]),2))
#            y.append(math.pow(float(coordinate[k][1])-float(coordinate[l][1]),2))
#            z.append(math.pow(float(coordinate[k][2])-float(coordinate[l][2]),2))
#        for k in range(len(atom)*len(atom)):
#          d.append(math.sqrt(x[k]+y[k]+z[k]))    
#        d = list(filter(None,d))
#       print(min(d))
#        m = min(d)
#-----------------------------------------------------------------------------------------------
#        n = len(atom)  
#        energy[0] = -289.19635 - float(energy[0]) / n
        if energy <= 100 and energy > -200 : # and m > 4.8 : 
          i = i + 1
          xsf_name = 'Si36-' + str(i) + '.xsf'
          fp = open(xsf_name, 'w')
          fp.write("# total energy = " + str(energy[0]) + " eV\n")
          fp.write("\n")
          fp.write("ATOMS\n")
          num = 0
          for j in range(len(atom)):
            fp.write(atom[j] + "      ")
            for l in range(3):
                fp.write("      " + str(float(coordinate[num][l])*0.52917721092))
            for l in range(3):
                fp.write("      " + str(float(force[num][l])*(27.211385/0.52917721092)))
            fp.write("\n")
            num += 1
          fp.close()
          coordinate = []
          force = []
          atom = []
          energy = []
          x = []
          y = []
          z = []
          d = []
        else:
          coordinate = []
          force = []
          atom = []
          energy = []
          x = []
          y = []
          z = []
          d = []
  return 0

if __name__ == '__main__':
    read_c_f()
