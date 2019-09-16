list1 = {'key1':'1','key2':'2','key3':'3'}
i = 0
for a in list1:
  list2 = dict(list1)
  for b, c in list1.items():
    print(list2)
    list2.pop(b)
    i = i + 1
    print(i)