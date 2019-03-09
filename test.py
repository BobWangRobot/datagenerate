import iotbx
from AEV import Rpart

rs_values = [ 0.900000,1.437500, 1.975000, 2.512500, 3.050000, 3.587500, 4.125000, 4.662500]
ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
angular_rs_values = [0.900000, 2.200000]
radial_cutoff = 5.2
angular_cutoff = 3.5
radial_nu = 8
radial_angular_nu = 8
angular_zeta = 8
Rj = [2.1, 2.2, 2.5]

a = Rpart('acid060.pdb')
print(a.Req(Rj))


