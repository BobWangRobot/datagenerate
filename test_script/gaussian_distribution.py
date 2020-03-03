import numpy as np
import matplotlib.pyplot as plt
import math
import random

def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

def plot_GD():
  Rs_list = [2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0]
  font1 = {'family': 'Times New Roman',
           'weight': 'normal',
           'size': 20,
           }
  for Rs in Rs_list:
    R = np.linspace(0, 15, 300)
    print(R)
    n = 4.0
    print(Rs)
    mR = np.exp(- n * ((R - Rs) ** 2))
    print(mR)
    color = randomcolor()
    plt.plot(R, mR, color, label='Rs=%r'%Rs)
    plt.legend(prop=font1, loc='best')
  plt.title(" gaussian distribution ", font1)
  plt.xlabel("R values", font1)
  plt.ylabel("Rs values", font1)
  plt.xticks(R[::20], (range(0,30)))
  plt.show()

def main():
  plot_GD()

if __name__ == '__main__':
  main()


