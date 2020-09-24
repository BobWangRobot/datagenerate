import os



def index_path(num):
  workdir = "/home/bob/Desktop/datagenerate/test6/per_structure/per%s"%num
  indexpath = "/home/bob/Desktop/datagenerate/test6/per_structure/per%s/INDEX"%num
  for root,dirs,files in os.walk(workdir):
    sorted(files)
    for file in files:
      f=open(indexpath,'a')
      filename=os.path.join(root,file)
      f.write(filename+"\n")
      f.close()
if __name__ == '__main__':
    for i in range(0, 51, 5):
      if i != 0:
        index_path(i/10)