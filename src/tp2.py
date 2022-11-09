import metnum as mt
import numpy as np

def leer_karateclub(path = "../karateclub_matriz.txt"):
  M = np.zeros(34*34).reshape(34, 34)
  with open(path,'r') as file:
    # reading each line 
    i = 0
    for line in file:
      line = line.split()
      for j in range(0, len(line)):

        M[i][j] = line[j]
      i+=1
  return M

if __name__ == '__main__':
    M = np.random.rand(5,5)
    #M = np.array([1,2,3,7,8,9,4,5,6]).reshape(3, 3)
    N = np.transpose(M)
    O = M + N
    #print(O)
    #M = mylib.add_one(M)
    lib_tp2.print_hey()
    num = O.shape[0]
    #print(lib_tp2.power_iteration(O, 1000, 0.000001))

    (eigvalNuestros, eigvecNuestros) = lib_tp2.deflation(O, num, 10000, 0.000000001)
    (eigvalNP, eigvecNP) = np.linalg.eigh(O)
    ordenNuestro = np.argsort(eigvalNuestros)
    ordenNP = np.argsort(eigvalNP)
    
    eigvalNuestrosO = eigvalNuestros[ordenNuestro]
    eigvecNuestrosO = eigvecNuestros[::, ordenNuestro]
    eigvalNPO = eigvalNP[ordenNP]
    eigvecNPO = eigvecNP[::, ordenNP]

    epsilon = 0.00000001
    print(np.allclose(eigvalNuestros[ordenNuestro], eigvalNP[ordenNP],epsilon))
    prodinter = np.zeros(num)
    for i in range(num):
        prodinter[i] = abs(np.dot(eigvecNuestrosO[::,i], eigvecNPO[::,i]))
    todosunos = np.ones(num)
    print(np.allclose(prodinter, todosunos, epsilon))

    #print(eigvalNuestrosO)
    #print(eigvecNuestrosO)

    #print(eigvalNPO)
    #print(eigvecNPO)

def leer_karateclub(path = "../karateclub_matriz.txt"):
  M = np.zeros(34*34).reshape(34, 34)
  with open(path,'r') as file:
    # reading each line 
    i = 0
    for line in file:
      line = line.split()
      for j in range(0, len(line)):

        M[i][j] = line[j]
      i+=1
  return M

if __name__ == '__main__':
    A = leer_karateclub("../datasets/karateclub_matriz.txt")
    a, e = mt.power_iteration(A.astype(float), 1000, 0.001)
    print(e)

