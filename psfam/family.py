'''
family.py
@author: Ben Reggio

The class family is used to hold information about commuting families of pauli strings. Families have three
representations. The first representation is a string representation: IY YI YY. The second is a list representation,
where a list p[i] represents the strings x_i z_p[i] or z_i x_p[i]. This would be [1,2,3]. The final representation
is a matrix represnetation, where the matrix acting on the binary representation of i gives p[i]. This would be:
A = [1,0]
    [0,1]
Since A v_1 = v_1; and Av_2 = v_2. This returns the correct p[i]. This family should have access to all three
representations. The list representation will be used in the constructor. This class should be able to get all
other representations from this one.

The other functionality of this class is that it should be able to diagonalize the family. A family with a Matrix
representation A has a corresponding A^{N/2 % N-1}. When this is entered using set_q(q), the family object
will be able to diagonalize itself. The family will take the form U D U dagger, where D represents a diagonal
matrix which is equal to a pauli string in the z family. In this case U = exp(sum(x_i)) where x_i will be the
x strings in the x_list parameter, which is filled when the set(q) method is called.

After diagonalization, this code integrates with qiskit with the method apply_to_circuit(circuit). It will
add onto the end of the circuit the clifford gate corresponding to the matrix U.
'''

import psfam.pauli_utils as util      #methods tobin and get_string are used from pauli_utils
from numpy import pi

#This class will contain the information for each family of commmuting Pauli strings. It should be able to give the 
#String representation of the family with to_string(), the list representation with get_permutation(), and the matrix
#representation using get_generating_matrix().
class family:
    
    #Constructor takes in the list representation of the family and the size of the strings.
    def __init__(self,p,m):
        self.p = p
        self.q = False
        self.m = m
        self.coefficients = "unassigned"
        self.signs = [0]*(2**m -1)
        
    #q is used to determine the post-rotation unitary matrix that diagonalizes this family. 
    #q should be the N/2 or 1/2 power of the generating matrix in permutation form. The generators of
    #this permutation are placed in x_list, so that U = exp(sum x_i) where x_i are the strings in x_list
    def set_q(self,q):
        self.q = q
        l=[]
        for i in range(self.m):
            l = l + [q[2**i-1]]
        self.x_list=l
        for i in range(1,2**self.m):
            self.get_sign(i)
        
    #Used to print out all of the strings in this family in string form
    def __str__(self):
        p = self.p
        ps = []
        m = self.m
        for i in range(len(p)):
            bj = util.tobin(p[i],m)
            bi = util.tobin(i+1,m)
            ps = ps + [util.get_string(bi,bj)]
        return ",".join(ps)

    #Used to print out all of the strings in this family in string form
    def to_string(self):
        p = self.p
        ps = []
        m = self.m
        for i in range(len(p)):
            bj = util.tobin(p[i],m)
            bi = util.tobin(i+1,m)
            ps = ps + [util.get_string(bi,bj)]
        return ps
    
    #Given a quantum circuit object circuit, this applies the correct rotation for U which diagonalizes the family
    def apply_to_circuit(self,circuit):
        apply_e_x_series(circuit,self.x_list,self.m)
    
    #Set the coefficients of each measurement
    def set_coefficients(self,c):
        self.coefficients = c
      
    #Get the coefficients of each measurement
    def get_coefficients(self):
        return self.coefficients
    
    def get_permutation(self):
        return self.p
    
    #Return the m x m symmetric matrix representing this family
    def get_generating_matrix(self):
        data = []
        for i in range(self.m):
            row = util.tobin(self.p[2**i - 1],self.m)
            data = [row] + data
        return data

    #Summary of this family
    def properties(self):
        s = f''
        s = s + 'Qubits: ' + str(self.m) + '\nPermutation:' + str(self.get_permutation())
        s = s + '\nCoefficients: ' + str(self.get_coefficients()) + '\nGenerating Matrix:'
        A = self.get_generating_matrix()
        for row in A:
            s = s + '\n' + str(row)
        return s

    #Print strings used to diagonalize family. The unitary U = exp(i pi/4 (sum P_i)) is used, where P_i are
    #the strings which result from this method
    def get_diagonalizing_strings(self):
        ds = []
        for i in range(self.m):
            j = 2**i - 1
            ds = ds + [util.get_string(util.tobin(0,self.m),util.tobin(self.q[j],self.m))]
        return ds

    def get_sign(self,i):
        pb = util.tobin(i,self.m)
        #print(self.to_string()[i-1])
        #print("S_ " + str(i) + " " + str(self.p[i-1]))
        totx = [0]*self.m
        s=0
        for j in self.x_list:
            jb = util.tobin(j,self.m)
            dtp = util.dotp(pb,jb)
            if(dtp%2 == 1):
                #totx = util.bin_add(totx,jb)
                s=s+1
        #print(util.dotp(pb,totx))
        #print(s)
        #print(totx)
        #print(util.tobin(self.p[i-1],self.m))
        if((util.dotp(pb,util.tobin(self.p[i-1],self.m)) + s) % 4 == 0):
            self.signs[i-1] = 1
        else:
            self.signs[i-1] = -1

#Extra classes for the x and z families. Should behave very similar to family, but won't break for not having
#a permutation that goes along with it.
class xfamily(family):
    def __init__(self,m):
        self.p = False
        self.q = False
        self.m = m
        self.coefficients = "unassigned"
        self.signs = [1]*(2**self.m - 1)

    def get_diagonalizing_strings(self):
        ds = []
        for i in range(self.m):
            s = ["I"]*self.m 
            s[i] = "Y"
            ds = ds + ["".join(s)]
        return ds
        
    def get_permutation(self):
        return False
    
    def apply_to_circuit(self,circuit):
        for i in range(self.m):
            circuit.u(-pi/2,0,0,i)
    
    def __str__(self):
        p = self.p
        ps = []
        m = self.m
        bi = util.tobin(0,m)
        for i in range(1,2**self.m):
            bj = util.tobin(i,m)
            ps = ps + [util.get_string(bi,bj)]
        return ",".join(ps)

    def to_string(self):
        p = self.p
        ps = []
        m = self.m
        bi = util.tobin(0,m)
        for i in range(1,2**self.m):
            bj = util.tobin(i,m)
            ps = ps + [util.get_string(bi,bj)]
        return ps
    
    def properties(self):
        s = f'x family'
        s = s + '\nQubits: ' + str(self.m)
        s = s + '\nCoefficients: ' + str(self.get_coefficients())
        return s
    
    def get_generating_matrix(self):
        return False
    
class zfamily(family):
    def __init__(self,m):
        self.p = False
        self.q = False
        self.m = m
        self.coefficients = "unassigned"
        self.signs = [1]*(2**self.m)
        
    def get_permutation(self):
        return False
    
    def apply_to_circuit(self,circuit):
        return circuit
    
    def __str__(self):
        p = self.p
        ps = []
        m = self.m
        bj = util.tobin(0,m)
        for i in range(0,2**self.m):
            bi = util.tobin(i,m)
            ps = ps + [util.get_string(bi,bj)]
        return ",".join(ps)

    def to_string(self):
        p = self.p
        ps = []
        m = self.m
        bj = util.tobin(0,m)
        for i in range(0,2**self.m):
            bi = util.tobin(i,m)
            ps = ps + [util.get_string(bi,bj)]
        return ps
    
    def properties(self):
        s = f'z family'
        s = s + '\nQubits: ' + str(self.m)
        s = s + '\nCoefficients: ' + str(self.get_coefficients())
        return s
    
    def get_generating_matrix(self):
        return False

    def get_diagonalizing_strings(self):
        return ["".join(["I"]*self.m)]



#Circuit represents a quantum_circuit object from qiskit. This method implements a rotation on this circuit.
#This rotation is the type exp(i pi/4  (Z_j)). This list l contains the qubits which are sigma_z. For example:
#if l = [1,4,5], then Z_j = IZIIZZ (since it has a Z in qubits 1 4 and 5 but not 0 2 and 3)
def apply_e_z(circuit,l):
    for i in range(1,len(l)):           #Start by applying a controlled not gate between the first z qubit
        circuit.cx(l[i],l[0])           #and every other qubit
    
    circuit.sdg(l[0])                     #Next, apply an Sdg gate to the first qubit
    
    for i in range(1,len(l)):           #Finally apply a second set of controlled not gates
        circuit.cx(l[i],l[0])

#Utility method for apply_e_x. The goal is to return a list of the qubits which need to be acted on.
def get_change_list(index,m):
    s = list(bin(index)[2:])        #get binary vector for the index
    l=[]
    first = False
    for i in range(len(s)):
        if(s[len(s)-i-1] == '1'):
            l = l + [i]             #If the ith qubit should be changed, add i to the list l
    return l

#Circuit represents a quantum_circuit object from qiskit. This method implements a rotation on this circuit.
#This rotation is the type exp(i pi/4  (X_index)) where index is an int, and X_index is the member of the X family
#which is represented by this vector. X_5 = XIX etc...
def apply_e_x(circuit,index,m): 
    l = get_change_list(index,m)        #Get the qubits which need to be acted on
    for i in l:
        circuit.u(pi/2,0,0,i)            #Apply y rotation to each of these qubits
    apply_e_z(circuit,l)                #Apply e^(i pi/4 Z_i) (This step involves entanglement)
    for i in l:
        circuit.u(-pi/2,0,0,i)             #Apply y rotation again

#Similar to the previous method apply_e_x, this method is used for repeated applications of this method.
#The rotation exp((i pi/4)(IX + XX)) should be U_y (exp(i pi/4) (IZ + ZZ)) U_y instead of the longer method
#U_y exp(pi/4)(IZ) U_y U_y (exp(pi/4)(ZZ)) U_y. This would be extremely inefficient
def apply_e_x_series(circuit,x_list,m):
    change_lists = []
    for x in x_list:
        change_lists = change_lists + [get_change_list(x,m)]  #Get qubits acted on
    
    #This step could be inefficient in other contexts. We Y rotate all qubits, since it is generally necessary for
    #this problem. However, we really should check which qubits need to be rotated if we want to reuse this method
    for i in range(m):
        circuit.u(-pi/2,0,0,i)                                    #U_y
    for l in change_lists:
        apply_e_z(circuit,l)                                    #exp(Z_j) for each j
    for i in range(m):
        circuit.u(pi/2,0,0,i)                                     #U_y