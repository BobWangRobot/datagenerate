import os

workdir="/home/bob/Desktop/datagenerate/test6/tst1"
indexpath="/home/bob/Desktop/datagenerate/test6/tst1/INDEX"


for root,dirs,files in os.walk(workdir):
  sorted(files)
  for file in files:
    f=open(indexpath,'a')
    filename=os.path.join(root,file)
    f.write(filename+"\n")
    f.close()
