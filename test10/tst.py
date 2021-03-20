import json
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
import os
import sys
import copy
from libtbx.utils import null_out
from scitbx.array_family import flex
import __init__ as aev

def run():
  with open("index", "r") as f:
    for line in f.readlines():
      data = json.loads(line)
      for key, value in data.items():
        result = []
        print(key)
        filename = key
        pdb_inp = iotbx.pdb.input(file_name=filename)
        model = mmtbx.model.manager(
          model_input=pdb_inp,
          log=null_out())
        a = aev.AEV(model=model)
        b = aev.compare(a)
        recs = aev.format_HELIX_records_from_AEV(b)
        for key1, value1 in value.items():
          for item in recs:
            min1 = 100
            min2 = 100
            recod_list = item.split()
            if key1 == recod_list[4]:
              for i in value1:
                num1 = abs(int(recod_list[5])-i[0])
                num2 = abs(int(recod_list[8])-i[1])
                if num1 < min1:
                  min1 = num1
                if num2 < min2:
                  min2 = num2
            if min1 < 6 and min2 < 6:
              fmt = "{0:>}  {1:>2}  {2:>2}  {3:>}  {4:>}  {5:>3}({6:>4}) {7:>}  {8:>}  {9:>3}" \
                    "({10:>4})  {11:>30} "
              result.append(fmt.format(recod_list[0],recod_list[1],recod_list[2],recod_list[3],recod_list[4],
                                       recod_list[5],min1,recod_list[6],recod_list[7],recod_list[8],min2,recod_list[9]))
            else:
              fmt = "{0:>}  {1:>2}  {2:>2}  {3:>}  {4:>}  {5:>3}({6:>4}) {7:>}  {8:>}  {9:>3}" \
                    "({10:>4})  {11:>30} "
              result.append(fmt.format(recod_list[0], recod_list[1], recod_list[2], recod_list[3], recod_list[4],
                                       recod_list[5], None, recod_list[6], recod_list[7], recod_list[8], None,
                                       recod_list[9]))

      print('\n'.join(result))
if __name__ == '__main__':
  run()



