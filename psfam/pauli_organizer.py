'''
Psfam.py
@author: Ben Reggio

The purpose of this code is to create an object Pauli_Organizer, which assists VQE by partitioning pauli strings into the minimum set of commuting families. The following things are required.

1 parameter m - representing the size of the pauli strings (number of characters in the string). This should be given by the user based on the size of the 

A symmetric matrix which is similar to a companion matrix. This is solved the class matrix_generator, which has a method get_generator_matrix(m) which takes m and returns the appropriate A. A^i gives the matrices which represent the different families

The object Pauli_Organizer contains family objects which have the information about these individual families.
They will be supplied with the following information from Pauli_organizer:
    - The size m of the strings
    - The list representation of the family
    - A representation of the correct diagonalization of the family through the set_q(q) method
    - The correct coefficients of each measurement if such a decomposition is fed to this object through the input_pauli_decomps() method

A user may call get_families() method to get a list of the 2**m + 1 family objects
'''


from psfam.matrix_generator import *        #To get A matrix which is symmetric and similar to a companion matrix
from psfam.family import *                                  #Contains the important family class

from psfam.pauli_utils import *      #Self-explanatory methods or methods which are not relevant to the main purpose of
                               #Pauli_organizer are stored in pauli_utils
    
#Main user focused class. This class will take an integer m as an input and generate a set of family objects
#The family objects will contain a perfect partition of the strings.
#
#A secondary function of this class is to assign the correct circuits to the families. While generating the families,
#this object uses the set_q function of the families to assign the correct post-circuit rotations to the family objects
#
#A third functionality is to take in a decomposition of a Hamiltonian into Pauli Strings, and output the correct
#coefficients of each measurement so that the VQE circuit can be evaluated in the minimum amout of steps
class PauliOrganizer:

    #Constructor initiates partitioning:
    def __init__(self,m):
        super().__init__()
        self.m = m
        self.u = False
        self.N = 2**m
        self.f = self.build_family_group()                              #Partition strings
        self.pauli_decomposition = self.initialize_decomposition()      #set pauli_decomposition to 0 so that it can
                                                                        #it can be set individually if necessary
    
        
    #Returns list of families
    def get_families(self):
        return self.f

    #Used for setting the entirety of a Hamilton decomposition into pauli_decomposition parameter
    def input_pauli_decomps(self,pauli_dict):
        coords = ([0]*self.m,[0]*self.m)                                    #coords gives the coordinates of cp
        cp = ['I']*self.m                                                   #cp (current pauli) begins as the identity
        self.pauli_decomposition[0][0] = pauli_dict["".join(cp)]            #inputs the coefficient of cp in the input pauli_dict
        for k in range(self.N ** 2):
            #print("coords = " + str(coords))
            #print('cp = ' + str("".join(cp)))
            coords = increment_coords(coords,self.m - 1)                    #Increase coords to next Z_2 value
            cp = increment_string(cp,self.m - 1)                            #Go to next string with next coords
            i = to_int(coords[0])
            j = to_int(coords[1])
            self.pauli_decomposition[i][j] = pauli_dict["".join(cp)]        #inputs next string decomposition into
                                                                            #pauli_decomposition
        
    #Use this to input the decomposition string by string
    def input_pauli_decomp(self,pauli_string,value):
        i,j = get_coords(pauli_string,self.m)
        self.pauli_decomposition[i][j] = value
    
    #Sets decomposition to 0 so that it can be entered individually
    def initialize_decomposition(self):
        pc = []
        for i in range(self.N):
            row = []
            for j in range(self.N):
                row = row + [0]
            pc = pc + [row]
        return pc
     
    #This function partitions the strings and creates the families.
    def build_family_group(self):
        A = get_generating_matrix(self.m)
        self.generating_matrix = A
        sol = process_A(A,self.m)
        f=[]
        for i in range(self.N - 1):
            fam = family(sol[i],self.m)
            q_i = 0
            if(i%2 == 1):
                q_i = (int)((i+1)/2) - 1
            else:
                q_i = (int)((self.N + i)/2) - 1
            fam.set_q(sol[q_i])
            f = f + [fam]
        f = f + [xfamily(self.m),zfamily(self.m)]
        return f
    
    #Print the properties of this solution
    def properties(self):
        s = f''
        s = s + 'Qubits: ' + str(self.m) +  '\nGenerating Matrix:'
        A = self.generating_matrix
        for row in A:
            s = s + '\n' + str(row)
        return s
    
    #Uses pauli_decomposition to calculate the coefficients of each measurement
    def calc_coefficients(self):
        f=self.f
        N = len(f) - 1
        for family_index in range(N - 1):
            p = f[family_index].p
            c=[]
            for state_index in range(0,N):
                s = 0
                for string_index in range(1,N):
                    v = eigenvalue(string_index,state_index)
                    effect = self.pauli_decomposition[string_index][p[string_index-1]]*v*f[family_index].signs[string_index-1]
                    s = s + effect
                c=c+[s]
            f[family_index].set_coefficients(c)
        x_index = N-1
        z_index = N
        c=[]
        for state_index in range(0,N):
            s = 0
            for string_index in range(1,N):
                v = eigenvalue(string_index,state_index)
                effect = self.pauli_decomposition[0][string_index]*v
                s = s + effect
            c=c+[s]
        f[x_index].set_coefficients(c)
        c=[]
        for state_index in range(0,N):
            s = 0
            for string_index in range(0,N):
                v = eigenvalue(string_index,state_index)
                effect = self.pauli_decomposition[string_index][0]*v
                s = s + effect
            c=c+[s]
        f[z_index].set_coefficients(c)
        return f

    def get_family_of_string(self,s):
        p = self.f[0].p
        a = 1
        u=[]
        if(self.u == False):
            for l in range(self.N -1):
                u = u + [a]
                a = p[a-1]
            self.u = u
        else:
            u = self.u
        i,j = util.get_coords(s)
        if(j == 0):
            return self.f[self.N]
        if(i == 0):
            return self.f[self.N -1]
        k = (u.index(j) - u.index(i) - 1) % (self.N -1)
        return self.f[k]

    def get_fam_index(self,s):
        p = self.f[0].p
        a = 1
        u=[]
        if(self.u == False):
            for l in range(self.N -1):
                u = u + [a]
                a = p[a-1]
            self.u = u
        else:
            u = self.u
        i,j = util.get_coords(s)
        if(j == 0):
            return self.N
        if(i == 0):
            return self.N - 1
        k = (u.index(j) - u.index(i) - 1) % (self.N -1)
        return k


#This function takes the generating matrix as an input, and outputs a list of each family in its permutation representation.
#For example, the one element in the list might be [2,3,1] representing the family XZ YX ZY.
def process_A(A,m):
  l1 = []                               #List representation of family 1
  for i in range(m):                    #Each ROW in A represents a generating string in the family
    t = A[m-i-1]                        #t represents the generating string
    l1 = l1 + [t]
    
    #Here, we use linearity of the solution and carefully order the elements in our list. l1[i] is the permutation P(i). Therefore, once we know l[i] and l[j], we know l[i + j]. Every time we add another element to the list, we can add it to each previous element.
    for j in range(len(l1)-1):
      l1 = l1 + [bin_add(l1[j],t)]
    
  for i in range(len(l1)):
    l1[i] = to_int(l1[i])                       #Reformat binary vectors as integers
  
  data = [l1]                                   #data will contain the entire output
  prev_p = l1.copy()                            #prev_p is the permutation of the previous family
  N = 2**m - 1
  for i in range(N-1):
    next_p = []                                 #next_p is the permutation of the next family (the one we are making)
    for j in range(N):
      next_p = next_p + [prev_p[l1[j]-1]]       #This statement is equivalent to P^(i+1)(x) = P(P^(i)(x)) 
    data = data + [next_p]
    prev_p = next_p.copy()
  return data

#Method to get the dot product of v_i and v_j which represent the vectors in z_2^m
#So eigenvalue(5,6) can be found by doing 5 -> (1,0,1) ; 6 -> (1, 1, 0). (1,0,1) dot (1,1,0) = 1 + 0 + 0 = 1 
def eigenvalue(i,j):
    ib = list(bin(i)[2:])
    jb = list(bin(j)[2:])
    a = len(ib)
    b = len(jb)
    m=max(a,b)
    ib = tobin(i,m)
    jb = tobin(j,m)
    s=0
    for k in range(m):
        if((ib[k] == 1) & (jb[k] == 1)):
            s = s + 1
    s = (s+1)%2
    if s == 0:
        s = -1
    return s