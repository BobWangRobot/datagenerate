def tst(result):
  sum = [1, 2, 3]
  value = [3, 4, 5]
  print("value1:", (sum, value, result))
  def ints(value1, sum, result):
    sum.append(4)
    value1.append(6)
    result = result + 1
    # print("value2:", (sum, value1, result))
  ints(value1=value, sum=sum, result=result)
  print("value2:", (sum, value, result))

if __name__ == '__main__':
    tst(3)