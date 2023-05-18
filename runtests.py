'''
runtests.py
@author: Andrew Lytle
Reorganized by Ben Reggio

Unit tests for psfam
'''

from psfam.pauli_organizer import PauliOrganizer
import functools,itertools
from scipy.stats import unitary_group
import scipy.sparse as sp
import qiskit
from qiskit.quantum_info.states import Statevector
from qiskit.circuit.random import random_circuit
import numpy as np
from numpy.random import normal
import operator

op_I = sp.eye(2)
op_Z = sp.dia_matrix([[1,0],[0,-1]])
op_X = sp.dia_matrix([[0,1],[1,0]])
op_Y = -1j * op_Z @ op_X
pauli_ops = { "I" : op_I, "Z" : op_Z, "X": op_X, "Y" : op_Y }

def sp_kron_dok(mat_A, mat_B): 
    """Kronecker (tensor) product of two sparse matrices.
    
    Args:
        mat_A, mat_B: 2d numpy arrays
    Returns:
        Sparse matrix in "dictionary of keys" format (to eliminate zeros)
    """
    return sp.kron(mat_A, mat_B, format = "dok")

def contribution_to_EV(m, fam, counts, nshots):
    cfs = fam.get_coefficients()
    #cfs = [_.real for _ in cfs]
    #print(f"{np.sum(np.array(cfs)).real = }")
    measurement = 0

    for i in range(2**m):
        bi = bin(i)[2:]
        while(len(bi)<m):
            bi = '0' + bi
        if(bi in counts):
            measurement += counts[bi]*cfs[i]/nshots
    return measurement

def to_pauli_vec(mat):
    """Pauli string decomposition of a matrix.

    Args:
        2d numpy array, mat
    Returns:
        Dictionary pauli_vec, e.g. {'XX': 0.5, 'XY': 0.5}
    """
    pauli_vec = {} # the dictionary we are saving

    mat_vec = np.array(mat).ravel()
    num_qubits = int(np.log2(np.sqrt(mat_vec.size)))

    for pauli_string in itertools.product(pauli_ops.keys(), repeat = num_qubits):
        # construct this pauli string as a matrix
        ops = [ pauli_ops[tag] for tag in pauli_string ]
        op = functools.reduce(sp_kron_dok, ops)

        # compute an inner product, same as tr(A @ B) but faster
        op_vec = op.transpose().reshape((1, 4**num_qubits))
        coefficient = (op_vec * mat_vec).sum() / 2**num_qubits
        #if coefficient != 0:
        pauli_vec["".join(pauli_string)] = coefficient

    return pauli_vec

def statedict_to_arr(statedict):
    "Convert state dictionary to numpy array."
    keys = list(statedict.keys())
    m = len(keys[0])
    #print(m)
    result = np.zeros(2**m, dtype='complex')
    bitstrings = itertools.product('01',repeat=m)
    for i, bitstring in enumerate(bitstrings):
        key = functools.reduce(operator.add, bitstring)
        if key in keys:
            result[i] = statedict[key]
        else:
            result[i] = 0.0
    return result

def direct_EV_eval(Hmat, state_circuit):
    psi = Statevector(state_circuit)
    state = statedict_to_arr(psi.to_dict())  # array.
    direct_eval = np.vdot(state,np.dot(Hmat,state))
    return direct_eval

def random_evals(N):
    "List of Gaussian distributed numbers [(0,1),...]."
    vals = []
    for i in range(N):
        vals.append(normal())
    return vals

def dot_all(Ms):
    "Dot product of [M1, M2, ...]"
    res = np.identity(Ms[0].shape[0])
    for M in Ms[::-1]:
        res = np.dot(M, res)
    return res

def hc(M):
    "Hermitian conjugate of M."
    return M.conj().T

def random_H(N):
    "Random NxN Hermitian matrix."
    evs = random_evals(N)
    D = np.diag(evs)
    U = unitary_group.rvs(N)
    H = dot_all([U, D, hc(U)])
    return H

def families_EV_eval(m, Hmat, state_circuit, backend, nshots):
    decomp = to_pauli_vec(Hmat)
    PO = PauliOrganizer(m)
    PO.input_pauli_decomps(decomp)
    PO.calc_coefficients()
    fams = PO.f
    total = 0
    for fam in fams:
        _circuit = state_circuit.copy()
        fam.apply_to_circuit(_circuit)
        job = backend.run(qiskit.transpile(_circuit, backend), shots=nshots)
        result = job.result()
        counts = result.get_counts(_circuit)
        #logging.info(f"{counts = }")
        contribution = contribution_to_EV(m, fam, counts, nshots).real
        #logging.info(f"{contribution = }")
        total += contribution
    return total

def random_unit_test(m):
    nshots = 1  # statevector_simulator backend.
    N = 2**m
    Hmat = random_H(N)
    nominal_depth = 4
    rc = random_circuit(m, nominal_depth, measure=False)
    direct_eval = direct_EV_eval(Hmat, rc)
    backend_sim = qiskit.Aer.get_backend('statevector_simulator')
    total = families_EV_eval(m, Hmat, rc, backend_sim, nshots)
    print(f"{total = }")
    print(f"{direct_eval = }")
    return np.isclose(total, direct_eval.real)

results = []
for m in range(1,6):  # Slows down for m >~ 5.
    for test in range(1, 101):
        print(f"{m = } , {test = }")
        passed = random_unit_test(m)
        results.append(passed)
        if passed:
            continue
        else:
            print(f"Test failed with {m = }")
            break
    if not passed:
        break
print(f"Performed {len(results)} tests, "
      f"all tests passed = {np.array(results).all()==True}")
