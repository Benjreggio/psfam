"""
matrix_generator.py
@author: Ben Reggio 

This code generates a m x m companion matrix and finds a symmetric matrix which is similar to that companion matrix in GF(2)

This is done in a few steps.

First, a companion matrix C is found with companion_matrix
Then, a symmetric matrix F is found with friend_matrix() such that FC^T = CF
M = N^T F N, such that M is symmetric and has all 1s on the diagonal with find_N()
Then, M = L L^T for some L, which is found using split()
Now N^T F N = L L^T so F = [(N^T)^-1 L ] [L^T N^-1] = R R^T where R = (N^T)^-1 L
Since F = R R^T, and F C = C^T F, R^-1 C R = A is symmetric and similar to C. R is our result
"""


import galois as gf     #Generates irreducible_poly for us
from psfam.pauli_utils import *
GF2 = gf.GF(2)          #Our vectors and matrices will be in GF(2)

#Generates a companion matrix from an irreducible polynomial of size m + 1. 
def companion_matrix(m):
  field = gf.GF(2**m)                               #initialize Galois field
  a = list(field.irreducible_poly.coeffs[0:m])       #Find irreducible polynomial
  #print(a)
  C=[]
  for i in range(m-1):                              #Populate companion matrix by row
    row = [GF2(0)]*m
    row[i+1] = GF2(1)
    C = C + [row]
  lastrow = a
  #for i in range(m):
  #  lastrow = lastrow + [a[m-i-1]]
  C = C + [lastrow]
  return C

#Generates a vector which is useful for generating friend matrix F:
def get_b(a):
  m=len(a)
  b = [GF2(1)]                              #b_1 = 1
  for i in range(1,m):
    bn = GF2(0)
    for k in range(0,i):
      bn = bn + a[m-i+k]*b[k]               #Formula for b_i = sum a_{m-i-k} b_k
    b = b + [bn]
    #print("b = " + str(b))
  return b

#Generates a symmetric matrix which satisfies FC^T = CF.
def friend_matrix(C):
  m = len(C)
  a = C[m-1]
  b=get_b(a)                        #Uses vector based off of the irreducible polynomial
  firstrow = [GF2(0)]*m
  firstrow[0] = GF2(1)
  F = []
  for i in range(0,m):
    nextrow = [GF2(0)]*m
    nextrow[m-1-i] = GF2(1)         #Second diagonal should be all ones
    for j in range(1,i+1):
      nextrow[m-1-i + j] =  b[j]
    F = F + [nextrow]
  return F

def find_N(f):
  m = len(f)
  flippers = []
  base = -1
  for i in range(m):
    if(f[i][i] == GF2(1)):
      if base == -1:
        base = i
    else:
      flippers = flippers + [i]
  N = []
  #print("Base = " + str(base))
  for i in range(m):
    nextrow = [GF2(0)]*m
    nextrow[i] = GF2(1)
    if(i == base):
      for j in flippers:
        nextrow[j] = GF2(1)
    N = N + [nextrow]
  return N

#For a symmetric matrix with 1s on the diagonal A, gives L where L L^T = A
def split(M):
  n = len(M)
  if(n == 1):
    return [[GF2(1)]]
  for i in range(n):
    M2 = []
    for j in range(n):
      if(j != i):
        nextrow = M[j].copy()
        nextrow.pop(i)
        M2 = M2 + [nextrow]

    if(is_invertible(M2)):

      L = split(M2)
      if(L):

        if(is_invertible(L)):
          L_prime = invert(L)

          b = M[i].copy()
          b.pop(i)

          eta = mtimes(L_prime,b)

          """
          This removed code is used to verify that a_0 = 1
          s = GF2(0)
          for j in range(len(eta)):
            s = s + eta[j]
          if(s == GF2(1)):
            print("Case: ")
            show_GF_mat(L)
            print()
            show_GF_mat(M)
            print()
            print(eta)
            print(i)
            show_GF_mat(M2)
            print()
            print()
          """

          LR = []
          for j in range(n):
            nextrow = []
            if (j<i):
              nextrow = L[j][0:i] + [GF2(0)] + L[j][i:n-1]
            if(j == i):
              nextrow = eta[0:i] + [GF2(1)] + eta[i:n-1]
            if (j>i):
              nextrow = L[j-1][0:i] + [GF2(0)] + L[j-1][i:n-1]
            LR = LR + [nextrow]
          return LR
        else:
            show_GF_mat(L)
            print("ERROR L NOT INVERTIBLE")
  return False


#Returns a symmetric matrix which is similar to a companion matrix of size m
def get_generating_matrix(m):
    C = companion_matrix(m)

    #print("C = ")
    #for row in C:
    #    print(row)

    F = friend_matrix(C)

    #print("F = ")
    #for row in F:
    #    print(row)

    N = find_N(F)

    #print("N = ")
    #for row in N:
    #    print(row)


    M = mtimes(transpose(N),mtimes(F,N))

    #print("M = ")
    #for row in M:
    #    print(row)

    L = split(M)

    #print("L = ")
    #for row in L:
    #    print(row)

    R=mtimes(transpose(invert(N)),L)
    #NF = mtimes(R,transpose(R))

    #print("NF = ")
    #for row in NF:
    #    print(row)

    A=mtimes(invert(R),mtimes(C,R))
    for i in range(m):
        for j in range(m):
            if(A[i][j] == GF2(1)):
                A[i][j] = 1
            else:
                A[i][j] = 0
    return A


#This method assists in testing. It checks whether a given matrix will work as a generating matrix for a family
def check(A):
  m = len(A)
  N = 2**m
  one = [0]*(m-1) + [1]
  v=one.copy()
  hits = [False]*N
  hits[to_int(v)-1] = True
  for i in range(N-2):
    r = mat_times_vec(A,v)
    r_int = to_int(r)
    if(hits[r_int-1]):
      return False
    else:
      hits[r_int-1] = True
    v = r.copy()
  r = mat_times_vec(A,v)
  return r == one


#TESTING
"""
C = companion_matrix(3)
print("C=")
for row in C:
    print(row)
    
F = friend_matrix(C)
print("F = ")

for row in F:
    print(row)
    
A=get_generating_matrix(5)

for row in A:
    print(row)
for i in range(3,20):
    M = get_generating_matrix(i)
    print(check(M))
"""