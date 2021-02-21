# Author: Pranav Tupe
# Reference: https://github.com/quantastica/qosf-mentorship


# Importing the neccessary libaries
import math
import numpy as np
from numpy.random import choice 

# Creating Exception for the simulator
class QExceptin(Exception):
    def __init__(self, arg):
        self.__err = arg

# Gates for the simulator

X = np.array([
[0, 1],
[1, 0]
])

H = np.array([
[1/np.sqrt(2), 1/np.sqrt(2)],
[1/np.sqrt(2), -1/np.sqrt(2)]
])

# Storing the gates in dictionary
gate_dictionary={"x":X,"h":H,"cx":X}

# Defining projection operator
P0x0 = np.array([
[1, 0],
[0, 0]
])

P1x1 = np.array([
[0, 0],
[0, 1]
])

# Simulator functions

# Function to generate list of combinations of output, taking number of bits as argument
def bit_list_generator(bits):
    if bits==1:
        return ["0","1"]
    else:
        bit_list = []
        for i in bit_list_generator(bits -1):
            bit_list = bit_list + [i+"0",i+"1"]
        return bit_list

# Initialize the ground state of number of given qubits
def get_ground_state(num_qubits):
    return [1]+[0]*(2**num_qubits -1)

# Get the combined state of qubits
def get_combined_state(*qn):
    combined_state = []
    if len(qn) == 1:
        return qn[0]
    for i in range(1,len(qn)):
        if i == 1:
            combined_state = np.kron(qn[0], qn[1])
        else:
            combined_state = np.kron(combined_state, qn[i])
    return combined_state

# Adding custom gate
def add_custom_gate(gate_name,gate_def_matrix):
    if type(gate_def_matrix) == type(X):
        gate_dictionary[gate_name]=gate_def_matrix
    else:
        raise QExceptin("Gate defining matrix should be of type numpy.ndarray")

# Getting the matrix operator of gate for the target qubit, depending upon total qubits and the gate to be applied
def get_operator(total_qubits, gate_unitary, target_qubits):
    gate_matrix = gate_dictionary[gate_unitary.lower()]
    operator= np.identity(2)
    operator1= np.identity(2)
    
    
    #If gate is for single qubit
    if len(target_qubits) == 1 :
        if target_qubits[0] == 1 and total_qubits == 1:
            return gate_matrix
        for i in range(1,total_qubits):
            if target_qubits[0] == 1:
                if i == 1:
                    operator = np.kron(gate_matrix,np.identity(2))
                else:
                    operator = np.kron(operator,np.identity(2))
            else:
                if i+1 == target_qubits[0]:
                 operator = np.kron(operator,gate_matrix)
                else:
                    operator = np.kron(operator,np.identity(2))
        return operator
    
    
    #If gate is for 2 qubits
    elif len(target_qubits) == 2 :
        for i in range(1,total_qubits):
            if target_qubits[0] == 1:
                if i == 1:
                    operator = np.kron(P0x0,np.identity(2))
                else:
                    operator = np.kron(operator,np.identity(2))
            else:
                if i+1 == target_qubits[0]:
                    operator = np.kron(operator,P0x0)
                else:
                    operator = np.kron(operator,np.identity(2))
        for i in range(1,total_qubits):
            if target_qubits[0] == 1:
                if i== 1 and i+1 == target_qubits[1]:
                    operator1 = np.kron(P1x1,gate_matrix)
                elif i == 1:
                    operator1 = np.kron(P1x1,np.identity(2))
                elif i+1 == target_qubits[1]:
                    operator1 = np.kron(operator1,gate_matrix)
                else:
                    operator1 = np.kron(operator1,np.identity(2))
            elif target_qubits[1] == 1:
                if i == 1 and i+1 == target_qubits[0]:
                    operator1 = np.kron(gate_matrix,P1x1)
                elif i == 1:
                    operator1 = np.kron(gate_matrix,np.identity(2))
                elif i+1 == target_qubits[0]:
                    operator1 = np.kron(operator1,P1x1)
                else:
                    operator1 = np.kron(operator1,np.identity(2))
            else:
                if i+1 == target_qubits[0]:
                    operator1 = np.kron(operator1,P1x1)
                elif i+1 == target_qubits[1]:
                    operator1 = np.kron(operator1,gate_matrix)
                else:
                    operator1 = np.kron(operator1,np.identity(2))
        return operator + operator1

# The core function of the simulator. It takes the initial state, applies matrix operator on it (depending on the specified gate and target qubit) and returns the final state of the circuit. 
def run_program(initial_state, program):
    if math.log2(len(initial_state)) ==  int(math.log2(len(initial_state))):
        total_qubits = int(math.log2(len(initial_state)))
        for p in program:
            initial_state = np.dot(get_operator(total_qubits,p["gate"],p["target"]),initial_state)
    else:
        raise QExceptin("State not initiated properly")
    return initial_state

# Choose element from state_vector using weighted random and return it's bit string
def measure_all(state_vector):
    probabilities = state_vector**2
    elements = bit_list_generator(int(math.log2(len(state_vector))))
    return choice(elements,1, p=probabilities)[0]



# Final Output
# Returns the final output in form of object which contains all the statistical information of bit string and how many times it occured. 
def get_counts(state_vector, num_shots):
    count = dict()
    for i in range(num_shots):
        state = measure_all(state_vector)
        if state in count:
            count[state] += 1
        else:
            count[state] = 1

    return count
# Note: We need to specify the number of shots to be fired 