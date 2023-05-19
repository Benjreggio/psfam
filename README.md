

The psfam package is a python package which organizes the complete set of pauli strings of any size into commuting families for VQE, finding a unitary matrix from clifford operators which diagonalizes all of the strings, and implementing this unitary matrix on a qiskit circuit.

psfam can be installed with 

```
pip install psfam
```

## Usage

Below are instructions on how to use the Psfam package. It will show the process of using a PauliOrganizer object to partition strings, then calculating an expectation value on a qiskit simulator.

You'll need numpy and qiskit, so lets import those:

```python
from psfam.pauli_organizer import *
import numpy as np
from qiskit import *
from qiskit.circuit.library import EfficientSU2
```

Pick the number of qubits, and assign this to $m$. Pauli_organizer takes only the number of qubits as an argument. In its constructor, it partitions all the families. We can look at the properties of the object using properties(). This will tell us what $A$ matrix was used to generate all of the families. 


```python
m=2
PO = PauliOrganizer(m)

print(PO.properties())

print()

for fam in PO.get_families():
    print(fam)
```

    Qubits: 2
    Generating Matrix:
    [1, 1]
    [1, 0]
    
    XZ,YX,ZY
    XY,ZX,YZ
    IY,YI,YY
    IX,XI,XX
    II,IZ,ZI,ZZ
    

The get\_families() function returns a list of all of the $2^m + 1$ family objects. The final 2 objects will be the x and z families. 

For normal families, there will be a permutation which generates the family, and the generating matrix for that family $A^i$ can be printed from that permutation. The following methods can be used to get the individual properties. Printing a family will show the strings in that family

The coefficients of the family will tell you how to calculate the measurement at the end. These would be the coefficients $\alpha_{i,j}$ of the contribution of this family to the expectation value of the Hamiltonian:

$$
\langle H \rangle = \sum_i \alpha_{i,j} |\langle \psi | U_i | \chi_j \rangle |^2
$$

Since we have not entered the decomposition of the Hamiltonian, these can't be calculated yet


```python
fam = PO.get_families()[0]

print(fam)
print(fam.get_permutation())
print(fam.get_generating_matrix())
print(fam.get_coefficients())
```

    XZ,YX,ZY
    [2, 3, 1]
    [[1, 0], [1, 1]]
    unassigned
    

This information can also be accessed using properties()


```python
print(fam.properties())
```

    Qubits: 2
    Permutation:[2, 3, 1]
    Coefficients: unassigned
    Generating Matrix:
    [1, 0]
    [1, 1]
    

This part can be skipped if you wish to calculate your coefficients differently. Suppose we have decomposed our Hamiltonian into pauli strings. This function I wrote to simulate this. It generates a dictionary with the decomposition of the hamiltonian into each Pauli string. $H = \sum \beta_i P_i$, and the dictionary is $\{ P_i : \beta_i \}$


```python
def gen_random_pauli_dict(m):
    pauli_dict = dict()
    cp = ['I']*m
    coords = ([0]*m,[0]*m)
    pauli_dict["".join(cp)] = np.random.random()
    for i in range(4 ** m):
        coords = increment_coords(coords,m - 1)
        cp = increment_string(cp,m - 1)
        pauli_dict["".join(cp)] = np.random.random()
    return pauli_dict

Hamiltonian_decomp = gen_random_pauli_dict(m)
```

Use the input_pauli_decomps() function to input this into the Hamiltonian. Then, the calc_coefficients() function calculates all the coefficients of the families, and returns the set of families for convenience. Now, when we read the properties of the family, the coefficients are calculated:


```python
PO.input_pauli_decomps(Hamiltonian_decomp)
f = PO.calc_coefficients()
fam1 = f[0]
print(fam1.properties())
```

    Qubits: 2
    Permutation:[2, 3, 1]
    Coefficients: [-1.1435602789371138, 0.7861666705463014, 0.5364746867002518, -0.1790810783094393]
    Generating Matrix:
    [1, 0]
    [1, 1]
    

Now lets make a sample measurement. 

The family object has a method apply_to_circuit() which automatically applies the correct unitary to the circuit. We will create a circuit which has a variational form, a unitary matrix, and a measurement at the end.


```python
qc = QuantumCircuit(m,m)

ansatz = EfficientSU2(m,su2_gates=['rx', 'y'], reps=1)
for parameter in ansatz.parameters:
    ansatz = ansatz.bind_parameters({parameter: np.random.random()*(2*np.pi)})
qc.compose(ansatz,inplace = True)
fam1.apply_to_circuit(qc)
qc.measure(range(m), range(m))

qc.draw('mpl')
```




    
![Cirucit diagram](https://github.com/Benjreggio/psfam/raw/main/output_15_0.png)
    



This code is straight from here: https://qiskit.org/documentation/tutorials/circuits/1_getting_started_with_qiskit.html. qiskit has methods to visualize the results.


```python
backend_sim = Aer.get_backend('qasm_simulator')
job_sim = backend_sim.run(transpile(qc, backend_sim), shots=1024)
result_sim = job_sim.result()
counts = result_sim.get_counts(qc)
print(counts)
qiskit.visualization.plot_histogram(counts)
```

    {'11': 933, '00': 44, '01': 24, '10': 23}
    




    
![Counts histogram](https://github.com/Benjreggio/psfam/raw/main/output_17_1.png)
    



Now that we have this result, we can use the counts of each result to get the probability of getting that result:

$$
P(010) = |\langle 010 | U | \psi \rangle |^2
$$

So we can just multiply these by the corresponding coefficients. The coefficients are in ascending binary order. The measurement for this family is calculated below:


```python
cfs = fam1.get_coefficients()
measurement = 0

for i in range(2**m):
    bi = bin(i)[2:]
    while(len(bi)<m):
        bi = '0' + bi
    if(bi in counts):
        measurement = measurement + counts[bi]*cfs[i]/1024

print(measurement)
```

    0.2523065901195194
    

Now this is repeated for each family and all of the results are added together
