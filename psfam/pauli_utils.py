"""
pauli_utils.py
@author Ben Reggio

This class contains methods which fall into miscallaneous categories, which are used by the classes:
matrix_generator.py
Psfam.py
family.py

These methods include methods which are self-explanatory, or methods distract from the main purpose of these classes
This code will show how these methods are made if one is looking for more detail.

The following methods are contained in this class:

transpose
invert
mtimes
is_invertible
get_coords
increment_coords
increment_string
get_progression
m_times
mat_times_vec
string_process
tobin
get_string
bin_add
to_int
enum_strings
"""


import galois as gf
GF2 = gf.GF(2)          #Some vectors and matrices will be in GF(2)


def dotp(a,b):
    m = len(a)
    s=0
    for i in range(m):
        s = s + a[i]*b[i]
    return s

#Util method to get transpose of a matrix
def transpose(M):
  M2 = []
  for i in range(len(M)):
    row = []
    for j in range(len(M)):
      row = row + [M[j][i]]
    M2 = M2 + [row]
  return M2

#Util method to invert a square GF(2) valued matrix
def invert(C):
  m = len(C)
  A=[]
  for i in range(m):                        #Construct augmented matrix
    row = [GF2(0)]*m
    row[i] = GF2(1)
    A = A + [C[i] + row]
  
  for i in range(m):                        #build upper triangular matrix
    j=i+1
    while(A[i][i] != GF2(1)):               #Swap rows such that the diagonal is 1
      if A[j][i] == GF2(1):
        t = A[j].copy()
        A[j] = A[i].copy()
        A[i] = t
      j=j+1

    for j in range(i+1,m):                  #Add rows below such that the bottom triangular area is 0
      if A[j][i] == GF2(1):
        for k in range(2*m):
          A[j][k] = A[i][k] + A[j][k]

  for i in range(m):                        #Clear upper triangular area to 0
    for j in range(i+1,m):
      if A[i][j] == 1:
        for k in range(2*m):
          A[i][k] = A[i][k] + A[j][k]
  
  B=[]
  for i in range(m):
    B = B + [A[i][m:]]                      #Select augmented side
  return B

#Util method to multiply matrices or vectors
def mtimes(M1,M2):
  if(type(M2[0]) == list):              #If matrix times matrix:
    M3 = []
    m = len(M1)
    for i in range(m):
      row = []
      for j in range(m):
        s=GF2(0)
        for k in range(m):              # sum M1_{i,k} M2_{k,j}
          s = s + M1[i][k]*M2[k][j]
        row = row + [s]
      M3 = M3 + [row]
    return M3
  else:                                 #If matrix times vector:
    r = []
    m = len(M1)
    for j in range(m):
      s=GF2(0)
      for k in range(m):
          s = s + M1[j][k]*M2[k]        #sum M1_{i,k} V_{k}
      r = r + [s]
    return r

#Util method to determine if a matrix is invertible. Goes through part of the process of inverting a matrix
def is_invertible(B):
  m = len(B)
  A = []
  for i in range(m):
    row = [GF2(0)]*m
    row[i] = GF2(1)
    A = A + [B[i].copy() + row]
  
  for i in range(m):
    j=i+1
    while(A[i][i] != GF2(1)):
      if j>=m:
        return False
      if A[j][i] == GF2(1):
        t = A[j].copy()
        A[j] = A[i].copy()
        A[i] = t
      j=j+1


    for l in range(i+1,m):
      if A[l][i] == GF2(1):
        for k in range(2*m):
          A[l][k] = A[i][k] + A[l][k]

  return True

#Utility method which gets the coordinates in the X-Z table of the pauli string. For example, XIY = XIX * IIZ = 5,1
def get_coords(pauli_input,m="unspecified"):
    ps = pauli_input
    if(type(ps) == int):
        ps=list(str(ps))
    if(type(ps) == str):
        ps = list(ps)
    i = []
    j = []
    if(type(ps) == list):
        if(m == "unspecified"):
            m = len(ps)
        else:
            if m != len(ps):
                print("ERROR: Wrong pauli size. m = " + str(m) + " the string " + str(pauli_input) + " has size " + str(len(ps)))
                return null
        ps = list(ps)

        if ps[0] in str_chars:
            for k in range(m):
                if(ps[k] == 'I') or (ps[k] == 'X'):
                    i = i + [0]
                else:
                    i = i + [1]
            for k in range(m):
                if(ps[k] == 'I') or (ps[k] == 'Z'):
                    j = j + [0]
                else:
                    j = j + [1]
        if ps[0] in int_chars:
            for k in range(m):
                if(ps[k] == '0') or (ps[k] == '2'):
                    i = i + [0]
                else:
                    i = i + [1]
            for k in range(m):
                if(ps[k] == '0') or (ps[k] == '1'):
                    j = j + [0]
                else:
                    j = j + [1]
    while(len(i)<m):
        i = [0] + i
    while(len(j)<m):
        j = [0] + j
    return to_int(i),to_int(j)

#Utility method to generate all strings from the coordinates of the X-Z table
tuple_progression = [(0,0),(1,0),(0,1),(1,1)]
def increment_coords(list_tuple,n):
    new_z,new_x = list_tuple
    ct = (new_z[n],new_x[n])
    for i in range(3):
        if ct == tuple_progression[i]:
            ct_z,ct_x = tuple_progression[i+1]
            new_x[n] = ct_x
            new_z[n] = ct_z
            return (new_z,new_x)
    ct_z,ct_x = tuple_progression[0]
    new_x[n] = ct_x
    new_z[n] = ct_z
    return increment_coords((new_z,new_x),n-1)

#Utility method to increment a string along a table. III -> IIZ -> IIX -> IIY -> IZI etc.
#This can be used to get every string in a specific order. n should be the size of the string.
def increment_string(s,n):
    new_s = s
    for i in range(3):
        #print("s = " + str(s))
        if s[n] == str_chars[i]:
            new_s[n] = str_chars[i+1]
            return new_s
    new_s[n] = "I"
    return increment_string(s,n-1)
    
    
int_chars = ['0','1','2','3']
#Utility method to multiply a matrix A with a vector v in Z2
def mat_times_vec(A,v):
    m = len(A)
    r = []
    for i in range(m):
        s=0
        for j in range(m):
            s = s + A[i][j] * v[j]
        r = r + [s%2]
    return r

str_chars=['I','Z','X','Y']
def string_process(s):
  ns = ''
  for i in range(len(s)):
    ns = ns + str_chars[(int)(s[i])]
  return ns

#Utility method to get the binary representation of an integer with length s.
def tobin(x,s):
    return [(x>>k)&1 for k in range(s-1,-1,-1)]

#Utility method to get the string given its coordinates on the X-Z table
def get_string(i,j):
  c = [str_chars[i[l] + 2*j[l]] for l in range(len(i))]
  return "".join(c)

#Utility method to get the bitwise sum of two integers (3 + 5 = 6)
def bin_add(a,b):
  c=[]
  for i in range(len(a)):
    if(a[i] == b[i]):
      c = c + [0]
    else:
      c = c + [1]
  return c

#Utility method for representing a binary vector as an integer
def to_int(a):
  m=len(a)
  s=0
  for j in range(m):
    if(a[m-j-1] == 1):
      s = s + 2**j
  return s

#Utility method to enumerate all of the pauli strings of length m.
def enum_strings(m):
    s1 = ['I']*m
    sl=[]
    for i in range(4**m):
        sl = sl + ["".join(s1)]
        increment_string(s1,m-1)
    return sl

#This is used for printing out a GF2 valued matrix in a readable manner
def show_GF_mat(mat):
    for row in mat:
        newrow = []
        for i in range(len(row)):
            if row[i] == GF2(0):
                newrow = newrow + [0]
            else:
                newrow = newrow + [1]
        print(newrow)