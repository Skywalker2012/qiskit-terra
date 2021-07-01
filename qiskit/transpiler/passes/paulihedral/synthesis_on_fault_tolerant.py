# ......
#
#
#


from qiskit.circuit.quantumcircuit import QuantumCircuit



def matching_operator_positions(pauli_str_x, pauli_str_y):
    positions = []
    for i,  x_op in enumerate(pauli_str_x[0]):
        if x_op == pauli_str_y[0][i] and x_op != "I":
            positions.append(i)
    return positions


def cancellation_potential(pauli_layers):
    layer_num = len(pauli_layers)
    cancellation_potential_array = [0]*layer_num
    for i in range(1, layer_num):
        positions = matching_operator_positions(pauli_layers[i][0][0], pauli_layers[i-1][0][0])
        cancellation_potential_array[i - 1] = len(positions)
    return cancellation_potential_array


def finding_max_matching_layer_pairs(cancellation_potential_array):
    layer_num = len(cancellation_potential_array)
    unmatched_layer_list = list(range(layer_num))
    pairs = []
    while True:
        max_val = -1
        max_id = 0
        for i in range(len(unmatched_layer_list) - 1):
            if unmatched_layer_list[i+1] == unmatched_layer_list[i] + 1:
                if cancellation_potential_array[unmatched_layer_list[i]] > max_val:
                    max_id = i
                    max_val = cancellation_potential_array[unmatched_layer_list[i]]
        if max_val == -1:
            for i in unmatched_layer_list:
                pairs.append((i,i))
            break
        else:
            pairs.append((unmatched_layer_list[max_id], unmatched_layer_list[max_id]+1))
            unmatched_layer_list = unmatched_layer_list[:max_id]+unmatched_layer_list[max_id+2:]
    
    pairs = sorted(pairs, key=lambda pairs : pairs[0])
    return pairs

def non_identity_positions(pauli_str):
    positions = []
    for i, op in enumerate(pauli_str[0]):
        if op != 'I':
            positions.append(i)
    return positions


def cnot_tree(non_id_positions, outer_qubits):
    if len(non_id_positions) == 0:
        return [], -1
    elif len(non_id_positions) == 1:
        return [], non_id_positions[0]
    
    inner_qubits = []
    for position in non_id_positions:
        if position not in outer_qubits:
            inner_qubits.append(position)
    if outer_qubits != []:
        r1 = []
        r2 = []
        for qubit in inner_qubits:
            if qubit > outer_qubits[-1]:
                r2.append(qubit)
            else:
                r1.append(qubit)
    else:
        r1 = []
        r2 = inner_qubits
    outercnot_qubits = outer_qubits + r2
    left_cnotset = []
    for i in range(len(outercnot_qubits)-1):
        left_cnotset.append((outercnot_qubits[i],outercnot_qubits[i+1]))
    root = outercnot_qubits[-1]
    if r1 != []:
        for i in range(len(r1)-1):
            left_cnotset.append((r1[i],r1[i+1]))
        left_cnotset.append((r1[-1], root))
    right_cnotset = left_cnotset[::-1]
    return right_cnotset, root


def syn_pauli_string(qc, qubit_num, string, right_cnotset, root, coeff):

    for k in range(qubit_num):
        if string[k] == 'X':
            qc.h(k)
        if string[k] == 'Y':
            qc.s(k)
            qc.h(k)
    for k in reversed(right_cnotset):
        qc.cx(k[0], k[1])
        
    if root >= 0:
        qc.rz(2*coeff, root)
        
    for k in right_cnotset:
        qc.cx(k[0], k[1])
    for k in range(qubit_num):
        if string[k] == 'X':
            qc.h(k)
        if string[k] == 'Y':                    
            qc.h(k)
            qc.sdg(k)

def break_layer(pauli_layers):
    broken_pauli_layers = []
    for pauli_layer in pauli_layers:
        max_pauli_block_size = 2
        for pauli_block in pauli_layer:
            if len(pauli_block) > max_pauli_block_size:
                max_pauli_block_size = len(pauli_block)
                
        if max_pauli_block_size == 2:    
            broken_pauli_layers.append(pauli_layer)
        else:
            temp_layer = [0] * (max_pauli_block_size-1)
            for i in range(max_pauli_block_size-1):
                temp_layer[i] = []
            
            for pauli_block in pauli_layer:
                for i, pauli_str in enumerate(pauli_block[:-1]):
                    temp_layer[i].append([pauli_str, pauli_block[-1]])
            broken_pauli_layers += temp_layer
    return broken_pauli_layers        


def opt_ft_backend(pauli_layers):
    qubit_num = len(pauli_layers[0][0][0][0])
    
    pauli_layers = break_layer(pauli_layers)
    
    qc = QuantumCircuit(qubit_num)
    
    cancellation_potential_array = cancellation_potential(pauli_layers)
    
    pairs = finding_max_matching_layer_pairs(cancellation_potential_array)
    
    for pair in pairs:
        if pair[0] == pair[1]:
            pauli_block = pauli_layers[pair[0]][0]
            pauli_str = pauli_block[0]
            
            non_id_positions = non_identity_positions(pauli_str)
            right_cnotset, root = cnot_tree(non_id_positions, [])
            syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
            
            for pauli_block in pauli_layers[pair[0]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
        else:
            matching_positions = matching_operator_positions(pauli_layers[pair[0]][0][0], pauli_layers[pair[1]][0][0])
            non_id_positions = non_identity_positions(pauli_layers[pair[0]][0][0])
            right_cnotset, root = cnot_tree(non_id_positions, matching_positions)
            syn_pauli_string(qc, qubit_num, pauli_layers[pair[0]][0][0][0], right_cnotset, root, pauli_layers[pair[0]][0][0][1]*pauli_layers[pair[0]][0][-1])
            
            for pauli_block in pauli_layers[pair[0]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
                    
            
            non_id_positions = non_identity_positions(pauli_layers[pair[1]][0][0])
            right_cnotset, root = cnot_tree(non_id_positions, matching_positions)
            syn_pauli_string(qc, qubit_num, pauli_layers[pair[1]][0][0][0], right_cnotset, root, pauli_layers[pair[1]][0][0][1]*pauli_layers[pair[1]][0][-1])
            
            for pauli_block in pauli_layers[pair[1]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])  

    return qc